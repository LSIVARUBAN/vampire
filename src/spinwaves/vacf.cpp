//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Lennon Sivaruban 2026. All rights reserved.
//
//   Email: lennon.sivaruban@postgrad.manchester.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <new>
#include <vector>

// VAMPIRE headers
#include "atoms.hpp"
#include "errors.hpp"
#include "material.hpp"
#include "sim.hpp"
#include "sld.hpp"
#include "spinwaves.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

// sw module headers
#include "internal.hpp"
#ifdef FFT
#include <fftw3.h>

namespace spinwaves{
namespace internal{

   // Fn to add products of each atom's velocity with its own velocity at a later time
   static void accumulate_vacf_block(const std::vector<double>& history,
                                    const int local_atoms, 
                                    const int capacity,
                                    const int frames, 
                                    const int overlap,
                                    const int max_lag, 
                                    const int fft_size,
                                    double* work, 
                                    fftw_complex* transform,
                                    fftw_plan forward, 
                                    fftw_plan backward,
                                    std::vector<double>& correlation){

      // store only the non-negative half of the Fourier coefficients, including zero and the Nyquist frequency (sampling frequency / 2)
      const int frequencies = fft_size/2 + 1;
      std::vector<double> power(3*frequencies, 0.0); // one sum of squared Fourier magnitudes for each direction

      for(int atom = 0; atom < local_atoms; atom++){
         for(int component = 0; component < 3; component++){

            // locate this atom's direction in the stored array  (frames entries in this row belong to the current block)
            const double* velocity = &history[(3*static_cast<size_t>(atom) + component)*capacity];
            std::fill(work, work + fft_size, 0.0); // fill unused part with zeroes
            std::copy(velocity, velocity + frames, work);

            fftw_execute(forward); // express the velocity sequence as frequency contributions
            for(int f = 0; f < frequencies; f++){
               // complex coefficient has real part [0] and imaginary part [1], calculate power as a sum of squares of Re and Im parts
               power[component*frequencies + f] += transform[f][0]*transform[f][0] + transform[f][1]*transform[f][1];
            }

            // the final frames will overlap with the start of the next block, so remove products whose two frames are in both to avoid double counting
            if(overlap > 0){
               std::fill(work, work + fft_size, 0.0); // fill non overalpping part with zeroes
               std::copy(velocity + frames - overlap, velocity + frames, work); // copy the tail of the current block to the start of the temp array
               fftw_execute(forward); // express tail as frequency contributions 
               for(int f = 0; f < frequencies; f++){
                  power[component*frequencies + f] -= transform[f][0]*transform[f][0] + transform[f][1]*transform[f][1]; // subtract tail contributions from power
               }
            }
         }
      }

      // inverse transform of a sum is the sum of the inverse transforms so transform accumulated powers once for each direction, instead of performing a separate inverse transform for every atom
      for(int component = 0; component < 3; component++){
         for(int f = 0; f < frequencies; f++){
            transform[f][0] = power[component*frequencies + f];
            transform[f][1] = 0.0; // Im part is 0 for squared magnitudes
         }
         fftw_execute(backward); // recover sums of velocity products at each frame separation

         for(int lag = 0; lag <= std::min(max_lag, frames - 1); lag++){ // furtherst pair is frames-1 intervals apart
            correlation[3*lag + component] += work[lag]/fft_size; // remove fft_size factor added by fftw
         }
      }
   }

} // end of internal namespace
} // end of spinwaves namespace
#endif

namespace spinwaves{

   // vacf function compares each atoms velocity with its own velocity at a later time
   // program first equilibrates the system, then records velocities during production and write the correlation versus time
   void run_phonon_vacf(){

      #ifndef FFT
         err::zexit("spinwaves:phonon-vacf requires FFTW");
      #else
         if(!sld::enabled || sim::integrator != sim::suzuki_trotter){
            err::zexit("spinwaves:phonon-vacf requires sim:integrator = spin-lattice");
         }

         const uint64_t sample_rate = internal::vacf_sample_rate;
         const double sample_time = mp::dt_SI*sample_rate; // time between recorded frames in seconds
         const uint64_t expected_frames = sim::total_time/sample_rate; 
         const double lag_steps = std::min(std::floor(internal::vacf_max_correlation_time/sample_time + 1.0e-9), static_cast<double>(expected_frames - 1)); // convert time separation into a number of recorded intervals
         const int max_lag = static_cast<int>(lag_steps);
         const int capacity = 2*max_lag + 1; // max_lag retained frames and max_lag+1 new frames

         int fft_size = 1;
         while(fft_size < 2*capacity - 1) fft_size *= 2; // power of two, 2N - 1 transform length 

         int local_atoms = atoms::num_atoms;
         #ifdef MPICF
            local_atoms = vmpi::num_core_atoms + vmpi::num_bdry_atoms;
         #endif
         const uint64_t number_of_atoms = vmpi::all_reduce_sum(static_cast<uint64_t>(local_atoms)); // total no physical atoms

         // store masses for com removal
         std::vector<double> masses(local_atoms);
         double total_mass = 0.0;
         for(int atom = 0; atom < local_atoms; atom++){
            masses[atom] = mp::material[atoms::type_array[atom]].mass;
            total_mass += masses[atom];
         }
         total_mass = vmpi::all_reduce_sum(total_mass); 

         const size_t channels = 3*static_cast<size_t>(local_atoms); // channel is the recorded velocity sequence for one atom in one direction, for each atom store xyz as three consecutive rows
         if(channels > std::numeric_limits<size_t>::max()/capacity/sizeof(double)){ // check the size calculation so allocation is sufficient
            err::zexit("VACF velocity buffer is too large");
         }
         const double memory_gb = channels*static_cast<double>(capacity)*sizeof(double)/1.0e9;
         zlog << zTs() << "VACF velocity buffer on rank " << vmpi::my_rank << ": " << memory_gb << " GB" << std::endl;
         std::vector<double> history;
         try{history.resize(channels*capacity, 0.0);}
         catch(const std::bad_alloc&){err::zexit("Unable to allocate VACF history");}
         
         std::vector<double> correlation(3*(max_lag + 1), 0.0);
         std::vector<double> work(fft_size, 0.0); // temp velocities before transformation, product sums afterwards
         fftw_complex* transform = fftw_alloc_complex(fft_size/2 + 1); // only half the complex coefficients are needed because the input velocities are real
         fftw_plan forward = fftw_plan_dft_r2c_1d(fft_size, work.data(), transform, FFTW_ESTIMATE); // convert real velocities to complex frequency coefficients
         fftw_plan backward = fftw_plan_dft_c2r_1d(fft_size, transform, work.data(), FFTW_ESTIMATE); // convert their squared magnitudes back to real correlation sums 

         const double production_temperature = sim::temperature; 
         sim::temperature = sim::Teq;
         vout::data();
         while(sim::time < sim::equilibration_time){
            sim::integrate(std::min(sim::partial_time, sim::equilibration_time - sim::time));
            stats::update();
            vout::data();
         }
         sim::temperature = production_temperature; // restore production temperature
         stats::reset(); // reset stats to remove equilibration samples from averages

         int buffered_frames = 0; // occupied time slots
         uint64_t sampled_frames = 0; // total number of distinct production frames
         uint64_t production_steps = 0; // steps relative to the start of production
         std::vector<double> com_velocity(3, 0.0);

         while(production_steps < sim::total_time){
            uint64_t steps = std::min(sample_rate - production_steps%sample_rate, sim::partial_time - production_steps%sim::partial_time); // steps until the next sampling 
            steps = std::min(steps, sim::total_time - production_steps); // stop at production end
            sim::integrate(steps);
            production_steps += steps;

            if(production_steps%sample_rate == 0){
               std::fill(com_velocity.begin(), com_velocity.end(), 0.0);
               if(internal::vacf_remove_com){ // fill com velocity if removal is required otherwise keep zeroes
                  for(int atom = 0; atom < local_atoms; atom++){
                     com_velocity[0] += masses[atom]*atoms::x_velo_array[atom];
                     com_velocity[1] += masses[atom]*atoms::y_velo_array[atom];
                     com_velocity[2] += masses[atom]*atoms::z_velo_array[atom];
                  }
                  vmpi::all_reduce_sum(com_velocity);
                  for(int component = 0; component < 3; component++) com_velocity[component] /= total_mass; // mv/m
               }

               for(int atom = 0; atom < local_atoms; atom++){
                  const size_t offset = 3*static_cast<size_t>(atom)*capacity + buffered_frames; // start at atom's x row, add capacity moves to its y row, add 2*capacity moves to its z row
                  history[offset] = atoms::x_velo_array[atom] - com_velocity[0]; // subtract from history, not from the velocity arrays used by the integrator
                  history[offset + capacity] = atoms::y_velo_array[atom] - com_velocity[1];
                  history[offset + 2*capacity] = atoms::z_velo_array[atom] - com_velocity[2];
               }
               // count this frame as sampled and buffered
               buffered_frames++;
               sampled_frames++;

               // process a full buffer of frames at capacity
               if(buffered_frames == capacity){
                  internal::accumulate_vacf_block(history, local_atoms, capacity, buffered_frames, max_lag, max_lag, fft_size, work.data(), transform, forward, backward, correlation); // accumulate products of each atom's velocity with its own velocity at a later time
                  // move the max_lag frames to the start of each atom's row, leave the rest of the row for new samples
                  for(size_t channel = 0; channel < channels; channel++){
                     double* row = &history[channel*capacity];
                     std::copy(row + capacity - max_lag, row + capacity, row);
                  }
                  buffered_frames = max_lag; // new samples are appended after the retained frames
               }
            }
            // update statistics and output as usual
            if(production_steps%sim::partial_time == 0 || production_steps == sim::total_time){
               stats::update();
               vout::data();
            }
         }

         // add all pairs in the final block
         internal::accumulate_vacf_block(history, local_atoms, capacity, buffered_frames, 0, max_lag, fft_size, work.data(), transform, forward, backward, correlation);
         fftw_destroy_plan(forward); 
         fftw_destroy_plan(backward);
         fftw_free(transform); // free memory allocated by fftw
         vmpi::all_reduce_sum(correlation); // sum contributions from all ranks
         if(vmpi::my_rank != 0) return;

         // write the correlation to file
         std::ofstream output("phonon_vacf.dat"); 
         output << std::setprecision(17);
         output << "# Cx, Cy, Cz and Ctotal in A^2/ps^2\n";
         output << "# lag_ps Cx Cy Cz Ctotal Cnormalised Norigins\n";
         const double zero_lag = (correlation[0] + correlation[1] + correlation[2])/(number_of_atoms*static_cast<double>(sampled_frames)); // zero separation value is just mean squared velocity whcih is used to scale the final curve
         for(int lag = 0; lag <= max_lag; lag++){
            const uint64_t origins = sampled_frames - lag; // number of starting frames for this lag separation
            const double normalisation = number_of_atoms*static_cast<double>(origins); // for a separation of N lag frames, each atom supplies samplied_frames-N pairs, so divide by the total number of pairs for the average
            const double cx = correlation[3*lag]/normalisation;
            const double cy = correlation[3*lag + 1]/normalisation;
            const double cz = correlation[3*lag + 2]/normalisation;
            const double total = cx + cy + cz; // velocity dot product
            const double normalised = zero_lag > 0.0 ? total/zero_lag : std::numeric_limits<double>::quiet_NaN(); // divide by zero lag value to normalise the correlation to 1 at zero separation
            output << lag*sample_time*1.0e12 << " " << cx << " " << cy << " " << cz << " " << total << " " << normalised << " " << origins << "\n";
         }
         output.close(); 
         std::cout << "Written phonon_vacf.dat from " << sampled_frames << " velocity frames" << std::endl;
      #endif
   }

} // end of spinwaves namespace
