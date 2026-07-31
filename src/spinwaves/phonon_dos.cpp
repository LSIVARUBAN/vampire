//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   Dev: Lennon Sivaruban, 2026.
//
//   Email: lennon.sivaruban@postgrad.manchester.ac.uk
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

// Spin-wave module headers
#include "internal.hpp"

#ifdef FFT
#include <fftw3.h>
#endif

namespace {

   int phonon_dos_local_atoms = 0;
   int phonon_dos_samples = 0;
   std::uint64_t phonon_dos_total_atoms = 0;

   // use single precision for compatibility with large systems and to reduce memory usage
   std::vector<float> phonon_dos_velocity_samples;

   std::size_t phonon_dos_index(const int atom,
                                const int component,
                                const int sample){
      return (static_cast<std::size_t>(atom) * 3u +
              static_cast<std::size_t>(component)) *
             static_cast<std::size_t>(phonon_dos_samples) +
             static_cast<std::size_t>(sample);
   }

} // namespace

namespace spinwaves {

   bool phonon_dos_enabled(){
      return internal::phonon_dos;
   }

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
            MPI_Allreduce(&local_atoms,
                          &phonon_dos_total_atoms,
                          1,
                          MPI_UINT64_T,
                          MPI_SUM,
                          MPI_COMM_WORLD);
         #endif

         if(phonon_dos_samples < 2 || phonon_dos_total_atoms == 0){
            err::zexit("Cannot initialise phonon DOS with less than two samples or zero atoms");
         }

         const std::size_t local_series = static_cast<std::size_t>(phonon_dos_local_atoms) * 3u;
         if(local_series >
            std::numeric_limits<std::size_t>::max() /
            static_cast<std::size_t>(phonon_dos_samples)){
            err::zexit("Error - Memory allocation failure.");
         }

         const std::size_t stored_values = local_series * static_cast<std::size_t>(phonon_dos_samples);

         try{
            phonon_dos_velocity_samples.assign(stored_values, 0.0f);
         }
         catch(const std::bad_alloc&){
            err::zexit("Insufficient memory.");
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

   void record_phonon_dos_sample(const int sample){

      if(sample < 0 || sample >= phonon_dos_samples){
         err::zexit("Phonon-DOS index is outside the time series");
      }

      for(int atom = 0; atom < phonon_dos_local_atoms; atom++){
         phonon_dos_velocity_samples[phonon_dos_index(atom, 0, sample)] = static_cast<float>(atoms::x_velo_array[atom]);
         phonon_dos_velocity_samples[phonon_dos_index(atom, 1, sample)] = static_cast<float>(atoms::y_velo_array[atom]);
         phonon_dos_velocity_samples[phonon_dos_index(atom, 2, sample)] = static_cast<float>(atoms::z_velo_array[atom]);
      }
   }

   void calculate_phonon_dos(){

      #ifndef FFT
         err::zexit("spinwaves:phonon-dos requires VAMPIRE to be compiled with FFTW");
      #else

         const int number_of_frequencies = phonon_dos_samples / 2 + 1;
         std::vector<double> local_spectrum(number_of_frequencies, 0.0);
         std::vector<double> global_spectrum(number_of_frequencies, 0.0);

         double* time_series = fftw_alloc_real(phonon_dos_samples);
         fftw_complex* frequency_series = fftw_alloc_complex(number_of_frequencies);
         if(time_series == NULL || frequency_series == NULL){
            if(time_series != NULL) fftw_free(time_series);
            if(frequency_series != NULL) fftw_free(frequency_series);
            err::zexit("FFTW could not allocate the phonon-DOS work arrays");
         }

         fftw_plan plan = fftw_plan_dft_r2c_1d(phonon_dos_samples,
                                               time_series,
                                               frequency_series,
                                               FFTW_ESTIMATE);
         if(plan == NULL){
            fftw_free(time_series);
            fftw_free(frequency_series);
            err::zexit("FFTW could not create the phonon-DOS transform plan");
         }

         for(int atom = 0; atom < phonon_dos_local_atoms; atom++){
            for(int component = 0; component < 3; component++){
               for(int sample = 0; sample < phonon_dos_samples; sample++){
                  time_series[sample] =
                     static_cast<double>(
                        phonon_dos_velocity_samples[
                           phonon_dos_index(atom, component, sample)]);
               }

               fftw_execute(plan);

               for(int frequency = 0;
                   frequency < number_of_frequencies;
                   frequency++){
                  local_spectrum[frequency] +=
                     std::hypot(frequency_series[frequency][0],
                                frequency_series[frequency][1]);
               }
            }
         }

         fftw_destroy_plan(plan);
         fftw_free(time_series);
         fftw_free(frequency_series);

         // The velocity history is no longer required once the local spectra
         // have been accumulated, so release it before the MPI reduction.
         std::vector<float>().swap(phonon_dos_velocity_samples);

         #ifdef MPICF
            MPI_Reduce(&local_spectrum[0],
                       &global_spectrum[0],
                       number_of_frequencies,
                       MPI_DOUBLE,
                       MPI_SUM,
                       0,
                       MPI_COMM_WORLD);
         #else
            global_spectrum.swap(local_spectrum);
         #endif

         if(vmpi::my_rank == 0){
            std::ofstream output("phonon_dos_raw.dat");
            if(!output.is_open()){
               err::zexit("Could not open phonon_dos_raw.dat for phonon-DOS output");
            }

            const double sample_period =
               static_cast<double>(sim::partial_time) * mp::dt_SI;
            const double normalisation =
               3.0 * static_cast<double>(phonon_dos_total_atoms);

            output << "# Direct per-atom velocity FFT phonon DOS\n";
            output << "# atoms " << phonon_dos_total_atoms << "\n";
            output << "# samples " << phonon_dos_samples << "\n";
            output << "# sample_period_s " << std::setprecision(17)
                   << sample_period << "\n";
            output << "# frequency_THz mean_absolute_velocity_fft\n";
            output << std::setprecision(17);

            for(int frequency = 0;
                frequency < number_of_frequencies;
                frequency++){
               const double frequency_thz =
                  static_cast<double>(frequency) /
                  (static_cast<double>(phonon_dos_samples) * sample_period) /
                  1.0e12;
               output << frequency_thz << " "
                      << global_spectrum[frequency] / normalisation << "\n";
            }
            output.close();

            std::cout << "Wrote direct per-atom phonon DOS to phonon_dos_raw.dat"
                      << std::endl;
            zlog << zTs()
                 << "Wrote direct per-atom phonon DOS to phonon_dos_raw.dat"
                 << std::endl;
         }

      #endif
   }

} // namespace spinwaves
