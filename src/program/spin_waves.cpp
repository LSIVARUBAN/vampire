//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Joel Hirst 2022. All rights reserved.
//
//   Email: j.r.hirst@shu.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
// Standard Libraries
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

//------------------------------------------------------------------------------
// Program to calculate a simple time series
//------------------------------------------------------------------------------
void spin_waves(){

	// check calling of routine if error checking is activated
	if(err::check==true) std::cout << "program::spin_waves has been called" << std::endl;

	double temp=sim::temperature;

   // Set equilibration temperature only if continue checkpoint not loaded
   if(sim::load_checkpoint_flag && sim::load_checkpoint_continue_flag){}
   else{
	   // Set equilibration temperature
	   sim::temperature=sim::Teq;
   }

	// Output data
	vout::data();

	// Equilibrate system
	while(sim::time<sim::equilibration_time){

		sim::integrate(sim::partial_time);

		// Calculate magnetisation statistics
		stats::update();

		// Output data
		vout::data();
	}

   // Set temperature and reset stats only if continue checkpoint not loaded
   if(sim::load_checkpoint_flag && sim::load_checkpoint_continue_flag){}
   else{

      // set simulation temperature
	   sim::temperature = temp;

      // Reset mean magnetisation counters
      stats::reset();

   }

	// Perform Time Series
	const bool calculate_phonon_dos = spinwaves::phonon_dos_enabled();
	if(calculate_phonon_dos) spinwaves::initialise_phonon_dos();

	while(sim::time<sim::equilibration_time+sim::total_time){

		// Integrate system
		sim::integrate(sim::partial_time);

		// calculate spinwaves JRH
		//const int step = (sim::time-sim::equilibration_time)/sim::partial_time-1;

		// std::cout << "calling spinwave function." << std::endl;
		const int sample = (sim::time-sim::equilibration_time)/sim::partial_time-1;
		if(calculate_phonon_dos){
			spinwaves::record_phonon_dos_sample(sample);
		}
		else{
			spinwaves::fft_in_space(atoms::x_coord_array,
											atoms::y_coord_array,
											atoms::z_coord_array,
											sample);
		}

		// Calculate magnetisation statistics
		stats::update();

		// Output data
		vout::data();

	}

	if(calculate_phonon_dos) spinwaves::calculate_phonon_dos();
	else spinwaves::fft_in_time();

}

}//end of namespace program
