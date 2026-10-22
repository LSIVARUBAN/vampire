//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans and Andrea Meo 2014-2018. All rights reserved.
//
//-----------------------------------------------------------------------------
//

// Standard Libraries
#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>

// Vampire Header files
#include "atoms.hpp"
#include "errors.hpp"
#include "material.hpp"
#include "program.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "vmath.hpp"
#include "vmpi.hpp"
#include "spinwaves.hpp" // JRH

namespace program{

namespace{

//------------------------------------------------------------------------------
// Applied field for each time series stage and the alternating production field
//------------------------------------------------------------------------------
namespace applied_field{

// production field from the input (stored so not overwritten by equilibration field)
bool production_field_stored = false;
double production_field = 0.0; // T
double production_direction[3] = {0.0, 0.0, 1.0};

void store_production_field(){
	if(production_field_stored) return;
	production_field = sim::H_applied;
	for(int i=0;i<3;i++) production_direction[i] = sim::H_vec[i];
	production_field_stored = true;
}

// field sign after a number of reversals
double sign(const uint64_t reversals){
	return reversals%2 == 0 ? 1.0 : -1.0;
}

// threshold for the next reversal, alternating sign
double target(const uint64_t reversals){
	return sign(reversals)*sim::alternating_field_trigger;
}

// check the trigger for the next reversal, the comparison swaps between leq and geq with each reversal
bool trigger_reached(const double value, const uint64_t reversals){
	if(std::isnan(value)) return false;
	const bool leq = (sim::alternating_field_trigger_method == sim::end_condition_leq) == (reversals%2 == 0);
	return leq ? value <= target(reversals) : value >= target(reversals);
}

// select the equilibration field if requested before equilibration_time,
// otherwise the production field with the sign left by the reversals already applied
void select_field(const bool ending_step){
	store_production_field();
	const bool equilibrating = sim::time < sim::equilibration_time || (ending_step && sim::time == sim::equilibration_time && sim::time > 0);
	sim::H_applied = production_field;
	for(int i=0;i<3;i++) sim::H_vec[i] = production_direction[i];
	if(equilibrating && sim::equilibration_applied_field_set){
		sim::H_applied = sim::Heq;
		if(sim::equilibration_applied_field_unit_vector_set){
			for(int i=0;i<3;i++) sim::H_vec[i] = sim::Heq_vec[i];
		}
	}
	if(!equilibrating && sim::alternating_field_enabled) sim::H_applied = sign(sim::alternating_field_state.reversals)*production_field;
}

// log the alternating field schedule on rank zero
void log_schedule(){
	if(vmpi::my_rank != 0) return;
	std::stringstream message;
	if(sim::alternating_field_mode == sim::alternating_field_interval_mode){
		message << "Alternating applied field: reversing every " << sim::alternating_field_interval_steps << " time steps from time step " << sim::equilibration_time;
	}
	else{
		const char* first_method = sim::alternating_field_trigger_method == sim::end_condition_leq ? "<=" : ">=";
		const char* second_method = sim::alternating_field_trigger_method == sim::end_condition_leq ? ">=" : "<=";
		const int component = static_cast<int>(sim::alternating_field_component);
		message << "Alternating applied field: reversing " << sim::alternating_field_buffer << " time steps after magnetisation_" << "xyz"[component] << " " << first_method << " " << target(0) << ", then " << second_method << " " << target(1) << ", alternately";
	}
	std::cout << message.str() << std::endl;
	zlog << zTs() << message.str() << std::endl;
}

// reverse the applied field before the step starting at sim::time and log the reversal, the physical field is H_applied*H_vec
void reverse(){

	sim::alternating_field_state_t& state = sim::alternating_field_state;
	const double old_field = sim::H_applied; // T
	state.reversals++;
	state.pending = false;
	sim::H_applied = sign(state.reversals)*state.base_field;

	if(vmpi::my_rank==0){
		const double real_time_ps = sim::time*mp::dt_SI*1.0e12; // start of the first step in the new field
		std::stringstream message;
		message << "Alternating field reversal " << state.reversals << " at time step " << sim::time << " (" << real_time_ps << " ps): B = ("
		        << old_field*sim::H_vec[0] << ", " << old_field*sim::H_vec[1] << ", " << old_field*sim::H_vec[2] << ") T -> ("
		        << sim::H_applied*sim::H_vec[0] << ", " << sim::H_applied*sim::H_vec[1] << ", " << sim::H_applied*sim::H_vec[2] << ") T";
		std::cout << message.str() << std::endl;
		zlog << zTs() << message.str() << std::endl;
	}

}

// apply any reversal due at the current step, then return the number of steps to the next sample, reversal or end of production
uint64_t reverse_if_due_and_get_batch_steps(const uint64_t end_time){

	sim::alternating_field_state_t& state = sim::alternating_field_state;
	const uint64_t production_time = sim::time - sim::equilibration_time; // production starts at equilibration_time exactly

	// fixed interval: reverse when another whole interval has passed since production started
	if(sim::alternating_field_mode == sim::alternating_field_interval_mode && production_time/sim::alternating_field_interval_steps > state.reversals) reverse();

	// trigger: reverse when the buffer after detection has passed
	if(sim::alternating_field_mode == sim::alternating_field_trigger_mode && state.pending && sim::time >= state.switch_time) reverse();

	uint64_t next_time = sim::equilibration_time + (production_time/sim::partial_time + 1)*sim::partial_time;
	next_time = std::min(next_time, end_time);
	if(sim::alternating_field_mode == sim::alternating_field_interval_mode) next_time = std::min(next_time, sim::equilibration_time + (state.reversals + 1)*sim::alternating_field_interval_steps);
	if(state.pending) next_time = std::min(next_time, state.switch_time);

	return next_time - sim::time;

}

// latch the trigger at a production sample, the reversal follows after the buffer and later samples do not postpone it
void detect_trigger(){

	sim::alternating_field_state_t& state = sim::alternating_field_state;
	if(sim::alternating_field_mode != sim::alternating_field_trigger_mode || state.pending || sim::time <= sim::equilibration_time) return;

	const int component = static_cast<int>(sim::alternating_field_component);
	const double value = stats::system_magnetization.get_magnetization()[component];
	if(!trigger_reached(value, state.reversals)) return;

	state.pending = true;
	state.detection_time = sim::time;
	state.switch_time = sim::time + sim::alternating_field_buffer;

	if(vmpi::my_rank==0){
		std::cout << "Alternating field trigger reached at time step " << sim::time << " (magnetisation_" << "xyz"[component] << " = " << value << "). Reversing field at time step " << state.switch_time << std::endl;
		zlog << zTs() << "Alternating field trigger reached at time step " << sim::time << " (magnetisation_" << "xyz"[component] << " = " << value << "). Reversing field at time step " << state.switch_time << std::endl;
	}

}

} // end of namespace applied_field

//------------------------------------------------------------------------------
// Early termination of the production run
//------------------------------------------------------------------------------
namespace end_conditions{

// if sim:end-condition-nan then stop simulation if all components of the system magnetisation are NaN
bool magnetisation_is_nan(const std::vector<double>& magnetisation){
	return std::isnan(magnetisation[0]) && std::isnan(magnetisation[1]) && std::isnan(magnetisation[2]) && std::isnan(magnetisation[3]);
}

void report_nan_magnetisation(){
	if(vmpi::my_rank==0){
		std::cout << "Simulation stopped at time step " << sim::time << " because the system magnetisation contains NaN" << std::endl;
		zlog << zTs() << "Simulation stopped at time step " << sim::time << " because the system magnetisation contains NaN" << std::endl;
	}
}

bool stop_for_nan_magnetisation(){
	if(!sim::end_condition_nan_enabled) return false;
	if(!magnetisation_is_nan(stats::system_magnetization.get_magnetization())) return false;
	report_nan_magnetisation();
	return true;
}

// latch the magnetisation end condition at a production sample and bring the requested end forward to the current time plus the
// buffer/ The latch and final step are saved in checkpoints so a continued run stops at the same step.
void latch_magnetisation_condition(const uint64_t production_end_time, uint64_t& requested_end_time){

	if(!sim::end_condition_enabled || sim::end_condition_triggered || sim::time >= production_end_time) return;

	const int component = static_cast<int>(sim::end_condition_component);
	const double component_value = stats::system_magnetization.get_magnetization()[component];
	const bool trigger_reached = sim::end_condition_trigger_method == sim::end_condition_leq ? component_value <= sim::end_condition_trigger : component_value >= sim::end_condition_trigger;
	if(!trigger_reached) return;

	sim::end_condition_triggered = true;
	const uint64_t remaining_time = production_end_time - sim::time;
	if(sim::end_condition_buffer < remaining_time) requested_end_time = sim::time + sim::end_condition_buffer;
	sim::end_condition_end_time = requested_end_time;

	if(vmpi::my_rank==0){
		std::cout << "Simulation end condition triggered at time step " << sim::time << " (magnetisation_" << "xyz"[component] << " = " << component_value << "). Continuing to final simulation step: " << requested_end_time << std::endl;
		zlog << zTs() << "Simulation end condition triggered at time step " << sim::time << " (magnetisation_" << "xyz"[component] << " = " << component_value << "). Continuing to final simulation step: " << requested_end_time << std::endl;
	}

}

} // end of namespace end_conditions

} // end of anonymous namespace

//------------------------------------------------------------------------------
// Select the applied field for the current time series stage.
//------------------------------------------------------------------------------
void set_time_series_field(const bool ending_step){
	applied_field::select_field(ending_step);
}

//------------------------------------------------------------------------------
// Program to calculate a simple time series
//------------------------------------------------------------------------------
void time_series(){

	// check calling of routine if error checking is activated
	if(err::check==true) std::cout << "program::time_series has been called" << std::endl;

	const double temp = sim::temperature;
	const bool continuing = sim::load_checkpoint_flag && sim::load_checkpoint_continue_flag;

	// keep the production field from the input, which the alternating field schedule reverses
	applied_field::store_production_field();
	if(sim::alternating_field_enabled){
		sim::alternating_field_state.base_field = applied_field::production_field;
		for(int i=0;i<3;i++) sim::alternating_field_state.base_direction[i] = applied_field::production_direction[i];
	}

	// set the equilibration temperature unless continuing from a checkpoint
	if(!continuing) sim::temperature = sim::Teq;

	// select the field for the current stage, the equilibration field or, when continuing in production, the field left by earlier reversals
	applied_field::select_field(true);

	// Output data
	vout::data();
	if(end_conditions::stop_for_nan_magnetisation()) return;

	// Equilibrate system
	while(sim::time<sim::equilibration_time){

		// The production field is restored after this loop exits, so if equilibration_time is not a multiple of partial_time
		// the last block would run past equilibration_time with the equilibration field still set. SLD switches its
		// damping/noise/thermostat parameters to production on the step where sim::time reaches equilibration_time, so when an
		// equilibration field is set, shorten the last block to end exactly at equilibration_time so the field switches on the same step.
		// An alternating field also needs production to start exactly at equilibration_time, where its schedule begins.
		const bool exact_equilibration_end = sim::equilibration_applied_field_set || sim::alternating_field_enabled;
		const uint64_t steps = exact_equilibration_end ? std::min(sim::partial_time, sim::equilibration_time - sim::time) : sim::partial_time;
		sim::integrate(steps);

		// Calculate magnetisation statistics
		stats::update();

		// Output data
		vout::data();

		// stop simulation if magnetisation is nan
		if(end_conditions::stop_for_nan_magnetisation()) return;
	}

	// set the production applied field from the input, with the sign left by any reversals when continuing
	applied_field::select_field(false);

	// set the production temperature and reset the mean magnetisation counters unless continuing from a checkpoint
	if(!continuing){
		sim::temperature = temp;
		stats::reset();
	}

	if(sim::alternating_field_enabled) applied_field::log_schedule();

	// Perform Time Series with optional early termination according to sim:end-condition. The number of total time steps are the maximum production length even when an end condition is enabled.
	// A latched end condition restored from a checkpoint keeps its final step.
	const uint64_t production_end_time = sim::equilibration_time + sim::total_time;
	uint64_t requested_end_time = production_end_time;
	if(sim::end_condition_enabled && sim::end_condition_triggered) requested_end_time = std::min(sim::end_condition_end_time, production_end_time);
	while(sim::time<requested_end_time){

		// Integrate system, with an alternating field apply any reversal due now and stop at the next sample, reversal or end of production
		const uint64_t steps = sim::alternating_field_enabled ? applied_field::reverse_if_due_and_get_batch_steps(requested_end_time) : sim::partial_time;
		sim::integrate(steps);

		// a batch ending at a reversal between output samples only advances the field schedule
		if(sim::alternating_field_enabled && (sim::time - sim::equilibration_time)%sim::partial_time != 0 && sim::time < requested_end_time) continue;

		// Calculate magnetisation statistics
		stats::update();

		// latch the alternating field trigger and the end condition before output, so a checkpoint written at this sample includes them
		if(sim::alternating_field_enabled) applied_field::detect_trigger();
		end_conditions::latch_magnetisation_condition(production_end_time, requested_end_time);

		// Output data
		vout::data();

		// apply the optional nan termination check
		if(end_conditions::stop_for_nan_magnetisation()) return;

	}

}

}//end of namespace program
