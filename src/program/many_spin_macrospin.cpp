//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Lennon Sivaruban and Richard F L Evans 2026. All rights reserved.
//
//------------------------------------------------------------------------------
//

// Standard libraries
#include <cmath>
#include <iostream>
#include <vector>

// Vampire Header files
#include "atoms.hpp"
#include "errors.hpp"
#include "material.hpp"
#include "program.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "vmpi.hpp"
#include "internal.hpp"

namespace program{

struct easy_axis{
	double x;
	double y;
	double z;
	int state;
};

// return the cubic state index or -1 if below geofencing threshold
static int aligned_state_cubic(const double mx, const double my, const double mz, const double threshold){
	const double inv_sqrt3 = 1.0 / std::sqrt(3.0);
	const easy_axis axes[8] = {
		{+inv_sqrt3, +inv_sqrt3, +inv_sqrt3, 0},
		{-inv_sqrt3, +inv_sqrt3, +inv_sqrt3, 1},
		{+inv_sqrt3, -inv_sqrt3, +inv_sqrt3, 2},
		{-inv_sqrt3, -inv_sqrt3, +inv_sqrt3, 3},
		{+inv_sqrt3, +inv_sqrt3, -inv_sqrt3, 4},
		{-inv_sqrt3, +inv_sqrt3, -inv_sqrt3, 5},
		{+inv_sqrt3, -inv_sqrt3, -inv_sqrt3, 6},
		{-inv_sqrt3, -inv_sqrt3, -inv_sqrt3, 7},
	};

	double best_dot = -1.0;
	int best_state = -1;
	for(const easy_axis& a : axes){
		const double dot = mx*a.x + my*a.y + mz*a.z;
		if(dot > best_dot){
			best_dot = dot;
			best_state = a.state;
		}
	}

	return (best_dot >= threshold) ? best_state : -1;
}

// return uniaxial state index or -1 if below threshold
static int aligned_state_uniaxial(const double mx, const double my, const double mz, const double threshold){
	const double ax = vout::uniaxial_axis_x; // uniaxial easy axis components
	const double ay = vout::uniaxial_axis_y;
	const double az = vout::uniaxial_axis_z;
	const easy_axis axes[2] = {
		{+ax, +ay, +az, 0},
		{-ax, -ay, -az, 1},
	};

	double best_dot = -1.0;
	int best_state = -1;
	for(const easy_axis& a : axes){
		const double dot = mx*a.x + my*a.y + mz*a.z;
		if(dot > best_dot){
			best_dot = dot;
			best_state = a.state;
		}
	}

	return (best_dot >= threshold) ? best_state : -1;
}

//------------------------------------------------------------------------------
// Program to calculate many spin macrospin switching statistics
//------------------------------------------------------------------------------
void ms_macrospin(){

	// check calling of routine if error checking is activated
	if(err::check==true) std::cout << "program::ms_macrospin has been called" << std::endl;

	double temp = sim::temperature;

	// Set equilibration temperature only if continue checkpoint not loaded
	if(sim::load_checkpoint_flag && sim::load_checkpoint_continue_flag){}
	else{
		sim::temperature = sim::Teq;
	}

	// Equilibrate system
	while(sim::time < sim::equilibration_time){

		sim::integrate(sim::partial_time);

        // Calculate magnetisation statistics
		stats::update();
	}

	// Set temperature and reset stats only if continue checkpoint not loaded
	if(sim::load_checkpoint_flag && sim::load_checkpoint_continue_flag){}
	else{

        // set simulation temperature
		sim::temperature = temp;

        // Reset mean magnetisation counters
		stats::reset();
	}

	const int num_atoms = atoms::num_atoms;

	const std::vector<double>& thresholds = vout::ms_macrospin_geofencing_thresholds;
	const size_t num_thresholds = thresholds.size();

	// track last aligned state per atom and count state transitions per threshold
	std::vector<std::vector<int>> last_aligned_state(num_thresholds, std::vector<int>(num_atoms, -1));
	std::vector<std::vector<uint64_t>> switch_counts(num_thresholds, std::vector<uint64_t>(num_atoms, 0u));
	std::vector<std::vector<uint64_t>> aligned_steps(num_thresholds, std::vector<uint64_t>(num_atoms, 0u));
	std::vector<std::vector<uint64_t>> lost_steps(num_thresholds, std::vector<uint64_t>(num_atoms, 0u));
	std::vector<std::vector<double>> first_transition_time(num_thresholds, std::vector<double>(num_atoms, 0.0));
	std::vector<std::vector<double>> last_transition_time(num_thresholds, std::vector<double>(num_atoms, 0.0));
	std::vector<std::vector<bool>> has_transition(num_thresholds, std::vector<bool>(num_atoms, false));

	const bool use_cubic = (program::internal::ms_macrospin_geofencing_mode == program::internal::ms_macrospin_geofencing_cubic);
	// Perform Time Series
	while(sim::time < sim::equilibration_time + sim::total_time){

        // Integrate system
		sim::integrate(sim::partial_time);

        // Calculate magnetisation statistics
		stats::update();

        // Geofence spins and count switches
		for(int atom = 0; atom < num_atoms; ++atom){

			if(!atoms::magnetic[atom]){
				continue;
			}

			const double mx = atoms::x_spin_array[atom];
			const double my = atoms::y_spin_array[atom];
			const double mz = atoms::z_spin_array[atom];

			for(size_t t = 0; t < num_thresholds; ++t){
				const double threshold = thresholds[t];
				const int aligned_state = use_cubic ? aligned_state_cubic(mx, my, mz, threshold) : aligned_state_uniaxial(mx, my, mz, threshold);

				// if aligned state is -1, the spin is lost
				if(aligned_state < 0){
					++lost_steps[t][atom];
					continue;
				}

				++aligned_steps[t][atom];

				// increase switch counts when aligned state changes
				if(last_aligned_state[t][atom] < 0){
					last_aligned_state[t][atom] = aligned_state;
				}
				else if(aligned_state != last_aligned_state[t][atom]){
					const double transition_time_s = sim::time * mp::dt_SI;
					if(!has_transition[t][atom]){
						first_transition_time[t][atom] = transition_time_s;
						has_transition[t][atom] = true;
					}
					last_transition_time[t][atom] = transition_time_s;
					++switch_counts[t][atom];
					last_aligned_state[t][atom] = aligned_state;
				}
			}
		}

		// calculate statistics, averaging tau over magnetic atoms, and output data at specified intervals
		vout::ms_macrospin_tau_avg.assign(num_thresholds, 0.0);
		vout::ms_macrospin_lost_time_avg.assign(num_thresholds, 0.0);
		vout::ms_macrospin_total_transitions.assign(num_thresholds, 0u);

		for(size_t t = 0; t < num_thresholds; ++t){
			uint64_t magnetic_atoms = 0u;
			double tau_sum = 0.0;
			double lost_time_sum = 0.0;
			uint64_t transitions_sum = 0u;

			for(int atom = 0; atom < num_atoms; ++atom){
				if(!atoms::magnetic[atom]){
					continue;
				}

				const uint64_t transitions = switch_counts[t][atom];
				const uint64_t total_steps = aligned_steps[t][atom] + lost_steps[t][atom];
				const double tau = (transitions > 0u && has_transition[t][atom])
					? ((last_transition_time[t][atom] - first_transition_time[t][atom]) / static_cast<double>(transitions))
					: 0.0;
				const double fractional_lost_time = (total_steps > 0u)
					? (static_cast<double>(lost_steps[t][atom]) / static_cast<double>(total_steps)) * 100.0
					: 0.0;

				tau_sum += tau;
				lost_time_sum += fractional_lost_time;
				transitions_sum += transitions;
				++magnetic_atoms;
			}

			if(magnetic_atoms > 0u){
				vout::ms_macrospin_tau_avg[t] = tau_sum / static_cast<double>(magnetic_atoms);
				vout::ms_macrospin_lost_time_avg[t] = lost_time_sum / static_cast<double>(magnetic_atoms);
			}
			vout::ms_macrospin_total_transitions[t] = transitions_sum;
		}

		vout::data();
	}
}

} // end of namespace program
