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

// C++ standard library headers
#include <cstddef>
#include <cmath>
#include <cstdint>
#include <vector>

// Vampire headers
#include "atoms.hpp"
#include "material.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

namespace stats{
namespace per_spin_geofencing{
namespace internal{

   struct alignment_result_t{
      double best_dot;
      int best_state;
   };

   mode_t geofencing_mode = disabled;
   int num_tracked_atoms = 0;
   uint64_t sampled_steps = 0u;

   std::vector<std::vector<int>> last_aligned_state;
   std::vector<std::vector<uint64_t>> switch_counts;
   std::vector<std::vector<uint64_t>> lost_steps;
   std::vector<std::vector<double>> first_transition_time;
   std::vector<std::vector<double>> last_transition_time;


   // Find the closest cubic axis for one spin
   alignment_result_t best_alignment_cubic(const double mx, const double my, const double mz){
      // When <111> is hard the six <100> directions are the easy axes
      if(vout::cubic_geofencing_111_hard){
         const double axes[6][3] = {
            {+1.0,  0.0,  0.0},
            {-1.0,  0.0,  0.0},
            { 0.0, +1.0,  0.0},
            { 0.0, -1.0,  0.0},
            { 0.0,  0.0, +1.0},
            { 0.0,  0.0, -1.0}
         };

         double best_dot = -1.0;
         int best_state = -1;
         for(int state = 0; state < 6; ++state){
            const double dot = mx*axes[state][0] + my*axes[state][1] + mz*axes[state][2];
            if(dot > best_dot){
               best_dot = dot;
               best_state = state;
            }
         }
         return {best_dot, best_state};
      }

      // Otherwise use the eight <111> easy axes
      const double inv_sqrt3 = 1.0 / std::sqrt(3.0);
      const double axes[8][3] = {
         {+inv_sqrt3, +inv_sqrt3, +inv_sqrt3},
         {-inv_sqrt3, +inv_sqrt3, +inv_sqrt3},
         {+inv_sqrt3, -inv_sqrt3, +inv_sqrt3},
         {-inv_sqrt3, -inv_sqrt3, +inv_sqrt3},
         {+inv_sqrt3, +inv_sqrt3, -inv_sqrt3},
         {-inv_sqrt3, +inv_sqrt3, -inv_sqrt3},
         {+inv_sqrt3, -inv_sqrt3, -inv_sqrt3},
         {-inv_sqrt3, -inv_sqrt3, -inv_sqrt3}
      };

      double best_dot = -1.0;
      int best_state = -1;
      for(int state = 0; state < 8; ++state){
         const double dot = mx*axes[state][0] + my*axes[state][1] + mz*axes[state][2];
         if(dot > best_dot){
            best_dot = dot;
            best_state = state;
         }
      }
      return {best_dot, best_state};
   }

   // Find whether one spin is closest to the positive or negative uniaxial axis
   alignment_result_t best_alignment_uniaxial(const double mx, const double my, const double mz){
      const double dot = mx*vout::uniaxial_axis_x + my*vout::uniaxial_axis_y + mz*vout::uniaxial_axis_z;
      return {std::fabs(dot), (dot >= 0.0) ? 0 : 1};
   }

} // namespace internal

// Select the geofencing mode used for per-spin statistics
void set_mode(const mode_t requested_mode){
   internal::geofencing_mode = requested_mode;
   return;
}

// Return the selected per-spin geofencing mode
mode_t get_mode(){
   return internal::geofencing_mode;
}

// Return whether per-spin geofencing is enabled
bool is_enabled(){
   return internal::geofencing_mode != disabled;
}

// Clear data and allocate one tracker per threshold and local atom
void reset(){
   if(!is_enabled()){
      return;
   }

   const std::size_t num_thresholds = vout::per_spin_geofencing_thresholds.size();
   internal::num_tracked_atoms = vmpi::num_local_atoms;
   internal::sampled_steps = 0u;

   internal::last_aligned_state.assign(num_thresholds, std::vector<int>(internal::num_tracked_atoms, -1));
   internal::switch_counts.assign(num_thresholds, std::vector<uint64_t>(internal::num_tracked_atoms, 0u));
   internal::lost_steps.assign(num_thresholds, std::vector<uint64_t>(internal::num_tracked_atoms, 0u));
   internal::first_transition_time.assign(num_thresholds, std::vector<double>(internal::num_tracked_atoms, 0.0));
   internal::last_transition_time.assign(num_thresholds, std::vector<double>(internal::num_tracked_atoms, 0.0));

   vout::per_spin_geofencing_tau_avg.assign(num_thresholds, 0.0);
   vout::per_spin_geofencing_lost_time_avg.assign(num_thresholds, 0.0);
   vout::per_spin_geofencing_total_transitions.assign(num_thresholds, 0u);
   return;
}

// Initialise the per-spin geofencing data structures
void initialize(){
   reset();
   return;
}

// Sample all locally owned magnetic spins and record lost states and transitions
void update(){
   if(!is_enabled()){
      return;
   }

   const std::size_t num_thresholds = vout::per_spin_geofencing_thresholds.size();
   if(internal::last_aligned_state.size() != num_thresholds || internal::num_tracked_atoms != vmpi::num_local_atoms){
      reset();
   }
   ++internal::sampled_steps;

   for(int atom = 0; atom < internal::num_tracked_atoms; ++atom){
      if(!atoms::magnetic[atom]){
         continue;
      }

      const internal::alignment_result_t alignment = (internal::geofencing_mode == cubic)
         ? internal::best_alignment_cubic(atoms::x_spin_array[atom], atoms::y_spin_array[atom], atoms::z_spin_array[atom])
         : internal::best_alignment_uniaxial(atoms::x_spin_array[atom], atoms::y_spin_array[atom], atoms::z_spin_array[atom]);

      for(std::size_t threshold_index = 0; threshold_index < num_thresholds; ++threshold_index){
         const int aligned_state = (alignment.best_dot >= vout::per_spin_geofencing_thresholds[threshold_index])
            ? alignment.best_state
            : -1;

         if(aligned_state < 0){
            ++internal::lost_steps[threshold_index][atom];
            continue;
         }

         int& previous_state = internal::last_aligned_state[threshold_index][atom];
         if(previous_state < 0){
            previous_state = aligned_state;
         }
         else if(aligned_state != previous_state){
            const double transition_time_s = sim::time * mp::dt_SI;
            uint64_t& transitions = internal::switch_counts[threshold_index][atom];
            if(transitions == 0u){
               internal::first_transition_time[threshold_index][atom] = transition_time_s;
            }
            internal::last_transition_time[threshold_index][atom] = transition_time_s;
            ++transitions;
            previous_state = aligned_state;
         }
      }
   }
   return;
}

// Calculate averages and reduce local per-spin results across all MPI processes
void synchronize(){
   if(!is_enabled()){
      return;
   }

   const std::size_t num_thresholds = vout::per_spin_geofencing_thresholds.size();
   std::vector<double> tau_sums(num_thresholds, 0.0);
   std::vector<double> lost_time_sums(num_thresholds, 0.0);
   std::vector<uint64_t> transition_sums(num_thresholds, 0u);
   uint64_t magnetic_atoms = 0u;

   for(int atom = 0; atom < internal::num_tracked_atoms; ++atom){
      if(!atoms::magnetic[atom]){
         continue;
      }
      ++magnetic_atoms;

      for(std::size_t threshold_index = 0; threshold_index < num_thresholds; ++threshold_index){
         const uint64_t transitions = internal::switch_counts[threshold_index][atom];
         const double tau = (transitions > 0u)
            ? ((internal::last_transition_time[threshold_index][atom] - internal::first_transition_time[threshold_index][atom]) /
               static_cast<double>(transitions))
            : 0.0;
         const double fractional_lost_time = (internal::sampled_steps > 0u)
            ? (static_cast<double>(internal::lost_steps[threshold_index][atom]) / static_cast<double>(internal::sampled_steps)) * 100.0
            : 0.0;

         tau_sums[threshold_index] += tau;
         lost_time_sums[threshold_index] += fractional_lost_time;
         transition_sums[threshold_index] += transitions;
      }
   }

   #ifdef MPICF
      if(num_thresholds > 0u){
         MPI_Allreduce(MPI_IN_PLACE, &tau_sums[0], static_cast<int>(num_thresholds), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, &lost_time_sums[0], static_cast<int>(num_thresholds), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, &transition_sums[0], static_cast<int>(num_thresholds), MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
      }
      MPI_Allreduce(MPI_IN_PLACE, &magnetic_atoms, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
   #endif

   vout::per_spin_geofencing_tau_avg.assign(num_thresholds, 0.0);
   vout::per_spin_geofencing_lost_time_avg.assign(num_thresholds, 0.0);
   vout::per_spin_geofencing_total_transitions = transition_sums;

   if(magnetic_atoms > 0u){
      const double inverse_magnetic_atoms = 1.0 / static_cast<double>(magnetic_atoms);
      for(std::size_t threshold_index = 0; threshold_index < num_thresholds; ++threshold_index){
         vout::per_spin_geofencing_tau_avg[threshold_index] = tau_sums[threshold_index] * inverse_magnetic_atoms;
         vout::per_spin_geofencing_lost_time_avg[threshold_index] = lost_time_sums[threshold_index] * inverse_magnetic_atoms;
      }
   }
   return;
}

} // namespace per_spin_geofencing
} // namespace stats
