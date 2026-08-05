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
#include <cmath>
#include <iostream>

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

// if sim:end-condition-nan then stop simulation if all components of the system magnetisation are NaN
bool magnetisation_is_nan(const std::vector<double>& magnetisation){

	return std::isnan(magnetisation[0]) &&
	       std::isnan(magnetisation[1]) &&
	       std::isnan(magnetisation[2]) &&
	       std::isnan(magnetisation[3]);

}

void report_nan_magnetisation(const std::vector<double>& magnetisation){

	if(vmpi::my_rank==0){
		std::cout << "Simulation stopped at time step " << sim::time << " because the system magnetisation contains NaN" << std::endl;
		zlog << zTs() << "Simulation stopped at time step " << sim::time << " because the system magnetisation contains NaN" << std::endl;
	}

}

bool stop_for_nan_magnetisation(){

	if(!sim::end_condition_nan_enabled) return false;

	const std::vector<double>& magnetisation = stats::system_magnetization.get_magnetization();
	if(!magnetisation_is_nan(magnetisation)) return false;

	report_nan_magnetisation(magnetisation);
	return true;

}

} // end of anonymous namespace

//------------------------------------------------------------------------------
// Program to calculate a simple time series
//------------------------------------------------------------------------------
void time_series(){

	// check calling of routine if error checking is activated
	if(err::check==true) std::cout << "program::time_series has been called" << std::endl;

	double temp=sim::temperature;
	double applied_field=sim::H_applied;

   // Set equilibration conditions only if continue checkpoint not loaded
   if(sim::load_checkpoint_flag && sim::load_checkpoint_continue_flag){}
   else{
	   // Set equilibration temperature
	   sim::temperature=sim::Teq;

	   // if requested, equilibrate in a different field strength
	   if(sim::equilibration_applied_field_set) sim::H_applied=sim::Heq;
   }

	// Output data
	vout::data();
	if(stop_for_nan_magnetisation()) return;

	// Equilibrate system
	while(sim::time<sim::equilibration_time){

		sim::integrate(sim::partial_time);

		// Calculate magnetisation statistics
		stats::update();

		// Output data
		vout::data();

		// stop simulation if magnetisation is nan
		if(stop_for_nan_magnetisation()) return;
	}

   // Set temperature and reset stats only if continue checkpoint not loaded
   if(sim::load_checkpoint_flag && sim::load_checkpoint_continue_flag){}
   else{

      // set production temperature and applied field
	   sim::temperature = temp;
	   sim::H_applied = applied_field;

      // Reset mean magnetisation counters
      stats::reset();

   }

	// Perform Time Series with optional early termination according to sim:end-condition. The number of total time steps are the maximum production length even when an end condition is enabled.
	const uint64_t production_end_time = sim::equilibration_time + sim::total_time;
	uint64_t requested_end_time = production_end_time;
	bool end_condition_triggered = false;
	while(sim::time<requested_end_time){

		// Integrate system
		sim::integrate(sim::partial_time);

		// Calculate magnetisation statistics
		stats::update();

		// Output data
		vout::data();

		// apply the optional nan termination check
		if(stop_for_nan_magnetisation()) return;

		if(sim::end_condition_enabled && !end_condition_triggered && sim::time < production_end_time){
			const std::vector<double>& magnetisation = stats::system_magnetization.get_magnetization();
			const int component = static_cast<int>(sim::end_condition_component);
			const double component_value = magnetisation[component];
			const bool trigger_reached = sim::end_condition_trigger_method == sim::end_condition_leq ? component_value <= sim::end_condition_trigger : component_value >= sim::end_condition_trigger; // check if the trigger condition is met

			if(trigger_reached){
				end_condition_triggered = true;
				const uint64_t remaining_time = production_end_time - sim::time;
				if(sim::end_condition_buffer < remaining_time)
					requested_end_time = sim::time + sim::end_condition_buffer; // set the end time to the current time plus the buffer

				if(vmpi::my_rank==0){
					const char* component_names[3] = {"magnetisation_x", "magnetisation_y", "magnetisation_z"};
					std::cout << "Simulation end condition triggered at time step " << sim::time << " (" << component_names[component] << " = " << component_value << "). Continuing to final simulation step: " << requested_end_time << std::endl;
					zlog << zTs() << "Simulation end condition triggered at time step " << sim::time << " (" << component_names[component] << " = " << component_value << "). Continuing to final simulation step: " << requested_end_time << std::endl; 
				}
			}
		}

	}

}

}//end of namespace program
