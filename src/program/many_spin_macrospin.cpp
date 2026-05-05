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

static void write_ms_macrospin_output(const int num_atoms,
										const std::vector<uint64_t>& switch_counts,
										const std::vector<uint64_t>& aligned_steps,
										const std::vector<uint64_t>& lost_steps,
										const std::vector<double>& first_transition_time,
										const std::vector<double>& last_transition_time,
										const std::vector<bool>& has_transition,
										const double total_time_s){
	if(vmpi::my_rank != 0){
		return;
	}

    // ignore file and screen output settings and write the outputs for many spin macrospin program
	if(zmag.is_open()){
		zmag.close();
	}
	vout::file_output_list.clear();
	vout::screen_output_list.clear();
	vout::output_rate = 1;
	vout::data();

    // set precision for outputs
	if(vout::custom_precision){
		zmag.precision(vout::precision);
		if(vout::fixed) zmag.setf(std::ios::fixed, std::ios::floatfield);
	}

    // loop through atoms and output switching statistics for each magnetic atom to file
	for(int atom = 0; atom < num_atoms; ++atom){
		
        if(!atoms::magnetic[atom]){
			continue;
		}

        // output file has columns: atom id, tau, fractional lost time, final simulation time, and total transitions. tau is defined as the simulation time between the final and initial transitions divided by the total number of transitions, and the fractional lost time represents the fraction of time when the geofencing threshold is not met.
		const uint64_t transitions = switch_counts[atom];
		const uint64_t total_steps = aligned_steps[atom] + lost_steps[atom];
		const double tau = (transitions > 0u && has_transition[atom])
			? ((last_transition_time[atom] - first_transition_time[atom]) / static_cast<double>(transitions))
			: 0.0;
		const double fractional_lost_time = (total_steps > 0u)
			? (static_cast<double>(lost_steps[atom]) / static_cast<double>(total_steps)) * 100.0
			: 0.0;
		zmag << atom << "\t" << tau << "\t" << fractional_lost_time << "\t" << total_time_s
			 << "\t" << transitions << std::endl;
	}
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

	// track last aligned state per atom and count state transitions
	std::vector<int> last_aligned_state(num_atoms, -1);
	std::vector<uint64_t> switch_counts(num_atoms, 0u);
	std::vector<uint64_t> aligned_steps(num_atoms, 0u);
	std::vector<uint64_t> lost_steps(num_atoms, 0u);
	std::vector<double> first_transition_time(num_atoms, 0.0);
	std::vector<double> last_transition_time(num_atoms, 0.0);
	std::vector<bool> has_transition(num_atoms, false);

	const bool use_cubic = (program::internal::ms_macrospin_geofencing_mode == program::internal::ms_macrospin_geofencing_cubic);
	const double threshold = use_cubic ? vout::cubic_geofencing_threshold : vout::uniaxial_geofencing_threshold;

	const uint64_t output_steps = program::internal::ms_macrospin_output_steps; // number of steps between periodic outputs, 0 disables periodic outputs

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

			const int aligned_state = use_cubic ? aligned_state_cubic(mx, my, mz, threshold) : aligned_state_uniaxial(mx, my, mz, threshold);

			// if aligned state is -1, the spin is lost
			if(aligned_state < 0){
				++lost_steps[atom];
				continue;
			}

			++aligned_steps[atom];

			// increase switch counts when aligned state changes
			if(last_aligned_state[atom] < 0){
				last_aligned_state[atom] = aligned_state;
			}
			else if(aligned_state != last_aligned_state[atom]){
				const double transition_time_s = sim::time * mp::dt_SI;
				if(!has_transition[atom]){
					first_transition_time[atom] = transition_time_s;
					has_transition[atom] = true;
				}
				last_transition_time[atom] = transition_time_s;
				++switch_counts[atom];
				last_aligned_state[atom] = aligned_state;
			}
		}

        // output stats periodically if enabled
		if(output_steps > 0u && (sim::time % output_steps == 0u)){
			const double current_time_s = sim::time * mp::dt_SI;
			write_ms_macrospin_output(num_atoms,
			                           switch_counts,
			                           aligned_steps,
			                           lost_steps,
			                           first_transition_time,
			                           last_transition_time,
			                           has_transition,
			                           current_time_s);
		}
	}

    // count total time in seconds
	const double total_time_s = sim::time * mp::dt_SI;

	// output results to file at the end of the simulation
	write_ms_macrospin_output(num_atoms,
	                           switch_counts,
	                           aligned_steps,
	                           lost_steps,
	                           first_transition_time,
	                           last_transition_time,
	                           has_transition,
	                           total_time_s);
}

} // end of namespace program
