//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Richard F L Evans and Rory Pond 2016. All rights reserved.
//
//   Email: richard.evans@york.ac.uk and rory.pond@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
// Headers
#include "vio.hpp"

// vio module headers
#include "internal.hpp"

// Global output filestreams
std::ofstream zinfo;
std::ofstream zlog;
std::ofstream zmag;
std::ofstream zgrain;
std::ofstream dp_fields;

namespace vout{

   std::string output_file_name;
   std::string zLogProgramName; /// Program Name
   std::string zLogHostName; /// Host Name
   bool        zLogInitialised=false; /// Initialised flag
   #ifdef WIN_COMPILE
   	int      zLogPid; /// Process ID
   #else
   	pid_t    zLogPid; /// Process ID
   #endif

   bool custom_precision = false; // enable user selectable precision for data output
   unsigned int precision = 6; // variable to control output precision (digits)
   int fw_size = 11;
   int fw_size_int = 11;
   bool fixed = false; // fixed precision output
   bool header_option = false; // output column headers on output file
   int max_header=14;
	// Namespace variable declarations
	std::vector<unsigned int> file_output_list(0);
	std::vector<unsigned int> screen_output_list(0);

	// Variables to control rate of data output to screen, output file and grain file
	int output_rate=1;
	int output_grain_rate=1;
	//int output_screen_rate=1; needs to be implemented

	bool gnuplot_array_format=false;

   namespace grain{

      // internal variables
      int output_rate = 1000; // rate of output compared to calculation

      // grain output list
      std::vector<grain::output_t> output_list(0);
   }

   double cubic_geofencing_threshold = 0.95; // default dot product threshold
   int cubic_geofencing_material_id = -1; // use global magnetisation by default
   bool cubic_geofencing_output_dot = true; // output best dot product column by default
   bool cubic_geofencing_111_hard = false; // use <111> easy axes by default
   double uniaxial_geofencing_threshold = 0.95; // default dot product threshold
   int uniaxial_geofencing_material_id = -1; // use global magnetisation by default
   bool uniaxial_geofencing_output_dot = true; // output best dot product column by default
   double uniaxial_axis_x = 0.0; // uniaxial easy axis components, default 001
   double uniaxial_axis_y = 0.0;
   double uniaxial_axis_z = 1.0;

   std::vector<double> per_spin_geofencing_thresholds = {0.95}; // default threshold
   std::vector<double> per_spin_geofencing_tau_avg(1, 0.0);
   std::vector<double> per_spin_geofencing_lost_time_avg(1, 0.0);
   std::vector<uint64_t> per_spin_geofencing_total_transitions(1, 0u);
}
