//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Lennon Sivaruban, 2026.
//
//   Email: lennon.sivaruban@manchester.ac.uk
//
//------------------------------------------------------------------------------
//
//   Direct phonon density of states calculation from per-atom velocity FFTs.
//   The methods follow S. Nikolov et al., "Dissociating the phononic, magnetic
//   and electronic contributions to thermal conductivity: a computational study
//   in alpha-iron", J. Mater. Sci. 57, 10535-10548 (2022), Eqs. (1)-(3).
//
//   Each owned atom's velocity v_{a,i}(t_n), a = x,y,z, is recorded every
//   sim:time-steps-increment steps after equilibration, giving N samples
//   separated by dt_s = time-steps-increment * dt.
//
//   Each series is transformed with a real-to-complex FFT
//
//      V_{a,i}(k) = sum_{n=0}^{N-1} v_{a,i}(t_n) exp(-2 pi i k n / N)
//
//   The magnitudes are averaged over the 3M velocity components of M atoms
//
//      psi(k) = 1/(3M) sum_{i=1}^{M} ( |V_{x,i}(k)| + |V_{y,i}(k)| + |V_{z,i}(k)| )
//
//   at frequency f_k = k / (N dt_s), giving a resolution of 1/(N dt_s) and a
//   maximum (Nyquist) frequency of 1/(2 dt_s). psi is written unnormalised to
//   phonon_dos_raw.dat.
//
//------------------------------------------------------------------------------

// C++ standard library headers
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <new>
#include <vector>

// Vampire headers
#include "atoms.hpp"
#include "errors.hpp"
#include "material.hpp"
#include "sim.hpp"
#include "spinwaves.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

// SW module headers
#include "internal.hpp"

#ifdef FFT
#include <fftw3.h>
#endif

namespace {

   int phonon_dos_local_atoms = 0;
   int phonon_dos_samples = 0;
   std::uint64_t phonon_dos_total_atoms = 0;

   // single precision to reduce memory usage
   std::vector<float> phonon_dos_velocity_samples;

   // index of a velocity sample, stored as a time series [atom][component][sample]
   std::size_t phonon_dos_index(const int atom,
                                const int component,
                                const int sample){
      return (static_cast<std::size_t>(atom) * 3u + static_cast<std::size_t>(component)) * static_cast<std::size_t>(phonon_dos_samples) + static_cast<std::size_t>(sample);
   }

}

namespace spinwaves {

   // true when spinwaves:phonon-dos replaces the usual spectrum calculation
   bool phonon_dos_enabled(){
      return internal::phonon_dos;
   }

   // allocate the velocity history for all owned atoms over the production run
   void initialise_phonon_dos(){

      #ifndef FFT
         err::zexit("spinwaves:phonon-dos requires VAMPIRE to be compiled with FFTW");
      #else

         phonon_dos_samples = static_cast<int>(sim::total_time / sim::partial_time);

         #ifdef MPICF
            phonon_dos_local_atoms = vmpi::num_core_atoms + vmpi::num_bdry_atoms;
         #else
            phonon_dos_local_atoms = atoms::num_atoms;
         #endif

         std::uint64_t local_atoms = static_cast<std::uint64_t>(phonon_dos_local_atoms);
         phonon_dos_total_atoms = local_atoms;

         #ifdef MPICF
            MPI_Allreduce(&local_atoms, &phonon_dos_total_atoms, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
         #endif

         const std::size_t local_series = static_cast<std::size_t>(phonon_dos_local_atoms) * 3u;
         if(local_series > std::numeric_limits<std::size_t>::max() / static_cast<std::size_t>(phonon_dos_samples)){
            err::zexit("Error: Memory allocation failure.");
         }

         const std::size_t stored_values = local_series * static_cast<std::size_t>(phonon_dos_samples);

         try{
            phonon_dos_velocity_samples.assign(stored_values, 0.0f);
         }
         catch(const std::bad_alloc&){
            err::zexit("Error: Insufficient memory.");
         }

         const double local_memory_mb = static_cast<double>(stored_values * sizeof(float)) / 1.0e6;
         const double total_memory_gb = static_cast<double>(phonon_dos_total_atoms) * 3.0 * static_cast<double>(phonon_dos_samples) * sizeof(float) / 1.0e9;

         std::cout << "Phonon DOS rank " << vmpi::my_rank << " recording " << phonon_dos_local_atoms << " owned atoms using " << local_memory_mb << " MB" << std::endl;
         if(vmpi::my_rank == 0){
            std::cout << "Phonon DOS recording " << phonon_dos_total_atoms << " atoms for " << phonon_dos_samples << " samples... total velocity storage is approx " << total_memory_gb << " GB" << std::endl;
            zlog << zTs() << "Phonon DOS recording " << phonon_dos_total_atoms << " atoms for " << phonon_dos_samples << " samples... total velocity storage is approx " << total_memory_gb << " GB" << std::endl;
         }

      #endif
   }

   // store the current velocities of owned atoms as time sample n
   void record_phonon_dos_sample(const int sample){

      if(sample < 0 || sample >= phonon_dos_samples){
         err::zexit("Error: Phonon-DOS index is outside the time series");
      }

      for(int atom = 0; atom < phonon_dos_local_atoms; atom++){
         phonon_dos_velocity_samples[phonon_dos_index(atom, 0, sample)] = static_cast<float>(atoms::x_velo_array[atom]);
         phonon_dos_velocity_samples[phonon_dos_index(atom, 1, sample)] = static_cast<float>(atoms::y_velo_array[atom]);
         phonon_dos_velocity_samples[phonon_dos_index(atom, 2, sample)] = static_cast<float>(atoms::z_velo_array[atom]);
      }
   }

   // transform each velocity series, accumulate psi(k) and write it from rank 0
   void calculate_phonon_dos(){

      #ifndef FFT
         err::zexit("spinwaves:phonon-dos requires VAMPIRE to be compiled with FFTW");
      #else

         // k = 0 ... N/2 for a real input series
         const int number_of_frequencies = phonon_dos_samples / 2 + 1;
         std::vector<double> local_spectrum(number_of_frequencies, 0.0);
         std::vector<double> global_spectrum(number_of_frequencies, 0.0);

         double* time_series = fftw_alloc_real(phonon_dos_samples);
         fftw_complex* frequency_series = fftw_alloc_complex(number_of_frequencies);
         if(time_series == NULL || frequency_series == NULL){
            if(time_series != NULL) fftw_free(time_series);
            if(frequency_series != NULL) fftw_free(frequency_series);
            err::zexit("FFTW could not allocate the phonon-DOS arrays");
         }

         fftw_plan plan = fftw_plan_dft_r2c_1d(phonon_dos_samples, time_series, frequency_series, FFTW_ESTIMATE);
         if(plan == NULL){
            fftw_free(time_series);
            fftw_free(frequency_series);
            err::zexit("FFTW could not create the phonon-DOS plan");
         }

         // sum |V_{a,i}(k)| over owned atoms and components
         for(int atom = 0; atom < phonon_dos_local_atoms; atom++){
            for(int component = 0; component < 3; component++){
               for(int sample = 0; sample < phonon_dos_samples; sample++){
                  time_series[sample] = static_cast<double>(phonon_dos_velocity_samples[phonon_dos_index(atom, component, sample)]);
               }

               fftw_execute(plan);

               for(int frequency = 0; frequency < number_of_frequencies; frequency++){
                  local_spectrum[frequency] += std::hypot(frequency_series[frequency][0], frequency_series[frequency][1]);
               }
            }
         }

         fftw_destroy_plan(plan);
         fftw_free(time_series);
         fftw_free(frequency_series);

         // velocity history is not required anymore once the local spectra have been accumulated, release it before MPI reduction.
         std::vector<float>().swap(phonon_dos_velocity_samples);

         #ifdef MPICF
            MPI_Reduce(&local_spectrum[0], &global_spectrum[0], number_of_frequencies, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
         #else
            global_spectrum.swap(local_spectrum);
         #endif

         if(vmpi::my_rank == 0){
            std::ofstream output("phonon_dos_raw.dat");
            if(!output.is_open()){
               err::zexit("Error: Could not open phonon_dos_raw.dat.");
            }

            // dt_s in seconds and the 1/(3M) prefactor of psi(k)
            const double sample_period = static_cast<double>(sim::partial_time) * mp::dt_SI;
            const double normalisation = 3.0 * static_cast<double>(phonon_dos_total_atoms);

            output << "# Direct per-atom velocity FFT phonon DOS\n";
            output << "# atoms " << phonon_dos_total_atoms << "\n";
            output << "# samples " << phonon_dos_samples << "\n";
            output << "# sample_period_s " << std::setprecision(17) << sample_period << "\n";
            output << "# frequency_THz mean_absolute_velocity_fft\n";
            output << std::setprecision(17);

            // f_k = k / (N dt_s) in THz
            for(int frequency = 0; frequency < number_of_frequencies; frequency++){
               const double frequency_thz = static_cast<double>(frequency) / (static_cast<double>(phonon_dos_samples) * sample_period) / 1.0e12;
               output << frequency_thz << " " << global_spectrum[frequency] / normalisation << "\n";
            }
            output.close();

            std::cout << "Wrote direct per-atom phonon DOS to phonon_dos_raw.dat" << std::endl;
            zlog << zTs() << "Wrote direct per-atom phonon DOS to phonon_dos_raw.dat" << std::endl;
         }

      #endif
   }

}
