//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Mara Strungaru 2022. All rights reserved.
//
//   Email: mara.strungaru@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <string>
#include <iostream>

// Vampire headers
#include "sld.hpp"
#include "errors.hpp"
#include "vio.hpp"

// sld module headers
#include "internal.hpp"

namespace sld{

   //---------------------------------------------------------------------------
   // Function to process input file parameters for sld module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line){

      // Check for valid key, if no match return false
      std::string prefix="spin-lattice";
      std::string prefix2="phonon"; //NEW
      if(key!=prefix && key!=prefix2) return false;
      //------------------------------------------------------------------------
      // If spin-transport parameter is requested, then enable module
      //------------------------------------------------------------------------
      sld::enabled = true;
      sld::internal::enabled=true;

      //----------------------------------
      // Now test for all valid options
      //----------------------------------

      std::string test = "linear-pump";
      if (word == test) {
         sld::internal::linear_pump_enabled = true; // This is a boolean flag
         return true; // We successfully matched the keyword
      }

      test = "frequency";
      if (word == test) {
         double val = atof(value.c_str());
         // Check for units in the 'unit' string provided by the parser
         if (unit == "THz") {
            val *= 1.0e12;
         } else if (unit == "GHz") {
            val *= 1.0e9;
         } // else, assume base units (Hz)
         sld::internal::phonon_frequency = val;
         return true;
      }

      test = "pulse-start-time";
      if (word == test) {
         double val = atof(value.c_str());
         // Check for units in the 'unit' string provided by the parser
         if (unit == "ps") {
            val *= 1.0e-12;
         } else if (unit == "fs") {
            val *= 1.0e-15;
         } // else, assume base units (s)
         sld::internal::phonon_pulse_start_time = val;
         return true;
      }

      test = "pulse-end-time";
      if (word == test) {
         double val = atof(value.c_str());
         // Check for units in the 'unit' string provided by the parser
         if (unit == "ps") {
            val *= 1.0e-12;
         } else if (unit == "fs") {
            val *= 1.0e-15;
         } // else, assume base units (s)
         sld::internal::phonon_pulse_end_time = val;
         return true;
      }

      // --- PARSING FOR FORCE AMPLITUDE ---
      test = "force-amplitude-x";
      if (word == test) {
         sld::internal::phonon_force_amplitude[0] = atof(value.c_str());
         return true;
      }
      test = "force-amplitude-y";
      if (word == test) {
         sld::internal::phonon_force_amplitude[1] = atof(value.c_str());
         return true;
      }
      test = "force-amplitude-z";
      if (word == test) {
         sld::internal::phonon_force_amplitude[2] = atof(value.c_str());
         return true;
      }

      // --- PARSING FOR WAVE LAMBDA ---
      test = "wave-lambda-x";
      if (word == test) {
         sld::internal::phonon_wave_lambda[0] = atof(value.c_str());
         return true;
      }
      //----------------------------------
      test = "wave-lambda-y";
      if (word == test) {
         sld::internal::phonon_wave_lambda[1] = atof(value.c_str());
         return true;
      }
      test = "wave-lambda-z";
      //----------------------------------
      if (word == test) {
         sld::internal::phonon_wave_lambda[2] = atof(value.c_str());
         return true;
      }
      // --- PARSING FOR WAVE DIRECTION ---
      test = "wave-direction-x";
      if (word == test) {
         sld::internal::phonon_wave_direction[0] = atof(value.c_str());
         return true;
      }
      //----------------------------------
      test = "wave-direction-y";
      if (word == test) {
         sld::internal::phonon_wave_direction[1] = atof(value.c_str());
         return true;
      }
      //----------------------------------
      test = "wave-direction-z";
      if (word == test) {
         sld::internal::phonon_wave_direction[2] = atof(value.c_str());
         return true;
      }
      //----------------------------------
      test = "potential";
      if( word == test ){
         test="harmonic";
         if( value == test ){
            sld::internal::lattice_potential = sld::internal::harmonic_lattice_potential;
            return true;
         }
         test="morse";
         if( value == test ){
            sld::internal::lattice_potential = sld::internal::morse_lattice_potential;
            return true;
         }
         test="snap";
         if( value == test ){
            sld::internal::lattice_potential = sld::internal::snap_lattice_potential;
            return true;
         }
         test="snapzbl";
         if( value == test ){
            sld::internal::lattice_potential = sld::internal::snap_zbl_lattice_potential;
            return true;
         }
      }

      // SNAP potentials have a coefficient and parameter file.
      test = "snap-coeff-file";
      if( word == test ){
         std::string snap_coeff_file=value;
         if(snap_coeff_file!=""){
            sld::internal::snap_potential.set_coeff_filename(snap_coeff_file);
            return true;
         }
         else{
            terminaltextcolor(RED);
            std::cerr << "Error - SNAP coefficient file not defined" << std::endl;
            terminaltextcolor(WHITE);
            return false;
         }
      }

      test = "snap-param-file";
      if( word == test ){
         std::string snap_param_file=value;
         if(snap_param_file!=""){
            sld::internal::snap_potential.set_param_filename(snap_param_file);
            return true;
         }
         else{
            terminaltextcolor(RED);
            std::cerr << "Error - SNAP parameter file not defined" << std::endl;
            terminaltextcolor(WHITE);
            return false;
         }
      }

      test = "snap-debug";
      if( word == test ){
         const bool debug = vin::check_for_valid_bool(value, word, line, prefix, "input");
         sld::internal::snap_potential.set_debug(debug);
         return true;
      }

      test = "zbl-inner-cutoff";
      if( word == test ){
         double r_c = vin::str_to_double(value);
         vin::check_for_valid_value(r_c, word, line, prefix, unit, "length", 0.1, 20.0,"input","0.1 - 20 A");
         sld::internal::zbl_inner_cutoff = r_c;
         return true;
      }

      test = "zbl-outer-cutoff";
      if( word == test ){
         double r_c = vin::str_to_double(value);
         vin::check_for_valid_value(r_c, word, line, prefix, unit, "length", 0.1, 20.0,"input","0.1 - 20 A");
         sld::internal::zbl_outer_cutoff = r_c;
         return true;
      }

      test = "zbl-atomic-number";
      if( word == test ){
         double z = vin::str_to_double(value);
         vin::check_for_valid_value(z, word, line, prefix, unit, "none", 1.0, 120.0,"input","1 - 120");
         sld::internal::zbl_atomic_number = z;
         return true;
      }

      test = "harmonic-debug";
      if( word == test ){
         const bool debug = vin::check_for_valid_bool(value, word, line, prefix, "input");
         sld::internal::harmonic_debug_enabled = debug;
         sld::internal::harmonic_debug_force_calls = 0;
         return true;
      }

      test = "coupling";
      if( word == test ){
         test="pseudodipolar";
         if( value == test ){
            sld::internal::pseudodipolar=true;
            sld::internal::full_neel=false;
            return true;
         }
         test="full-neel";
         if( value == test ){
            sld::internal::full_neel=true;
            sld::internal::pseudodipolar=false;
            return true;
         }
      }

      // inverse fourth or bethe-slater radial functions for the full Neel coupling
      test = "neel-radial-function";
      if( word == test ){
         if(value == "inverse-fourth"){
            sld::internal::neel_radial_function = sld::internal::inverse_fourth_neel_radial_function;
            return true;
         }
         if(value == "bethe-slater"){
            sld::internal::neel_radial_function = sld::internal::bethe_slater_neel_radial_function;
            return true;
         }
         err::zexit("spin-lattice:neel-radial-function must be inverse-fourth or bethe-slater");
      }

      // hard or smooth cutoff fns for the full neel coupling
      test = "neel-cutoff-function";
      if( word == test ){
         if(value == "hard"){
            sld::internal::neel_cutoff_function = sld::internal::hard_neel_cutoff_function;
            return true;
         }
         if(value == "smooth"){
            sld::internal::neel_cutoff_function = sld::internal::smooth_neel_cutoff_function;
            return true;
         }
         err::zexit("spin-lattice:neel-cutoff-function must be hard or smooth");
      }

      // Separate cutoff and switch distances for the full Neel coupling
      test = "neel-l-cutoff-range";
      if( word == test ){
         double cutoff = vin::str_to_double(value);
         vin::check_for_valid_value(cutoff, word, line, prefix, unit, "length", 0.1, 20.0, "input", "0.1 - 20 A");
         sld::internal::r_cut_neel_l = cutoff;
         sld::internal::r_cut_neel_l_set = true;
         return true;
      }

      test = "neel-q-cutoff-range";
      if( word == test ){
         double cutoff = vin::str_to_double(value);
         vin::check_for_valid_value(cutoff, word, line, prefix, unit, "length", 0.1, 20.0, "input", "0.1 - 20 A");
         sld::internal::r_cut_neel_q = cutoff;
         sld::internal::r_cut_neel_q_set = true;
         return true;
      }

      test = "neel-l-switch-range";
      if( word == test ){
         double switch_distance = vin::str_to_double(value);
         vin::check_for_valid_value(switch_distance, word, line, prefix, unit, "length", 0.1, 20.0, "input", "0.1 - 20 A");
         sld::internal::r_switch_neel_l = switch_distance;
         sld::internal::r_switch_neel_l_set = true;
         return true;
      }

      test = "neel-q-switch-range";
      if( word == test ){
         double switch_distance = vin::str_to_double(value);
         vin::check_for_valid_value(switch_distance, word, line, prefix, unit, "length", 0.1, 20.0, "input", "0.1 - 20 A");
         sld::internal::r_switch_neel_q = switch_distance;
         sld::internal::r_switch_neel_q_set = true;
         return true;
      }

      test = "exchange-function";
      if( word == test ){
         if(value == "cubic"){
            sld::internal::exchange_function = sld::internal::cubic_exchange_function;
            return true;
         }
         if(value == "bethe-slater"){
            sld::internal::exchange_function = sld::internal::bethe_slater_exchange_function;
            return true;
         }
         terminaltextcolor(RED);
         std::cerr << "Error - spin-lattice:exchange-function must be cubic or bethe-slater" << std::endl;
         terminaltextcolor(WHITE);
         return false;
      }

      test = "spin-hamiltonian";
      if( word == test ){
         if(value == "bilinear"){
            sld::internal::spin_hamiltonian =
               sld::internal::bilinear_spin_hamiltonian;
            return true;
         }
         if(value == "biquadratic"){
            sld::internal::spin_hamiltonian =
               sld::internal::biquadratic_spin_hamiltonian;
            return true;
         }
         terminaltextcolor(RED);
         std::cerr << "Error - spin-lattice:spin-hamiltonian must be bilinear or biquadratic" << std::endl;
         terminaltextcolor(WHITE);
         return false;
      }

      test = "exchange-offset";
      if( word == test ){
         // choose whether to shift the hamiltonian by -1
         sld::internal::exchange_offset =
            vin::check_for_valid_bool(value, word, line, prefix, "input");
         return true;
      }

      test = "spin-temperature-correction";
      if( word == test ){
         sld::internal::spin_temperature_correction =
            vin::check_for_valid_bool(value, word, line, prefix, "input");
         return true;
      }

      test = "thermostat"; // SLED or standard thermostat, if not specified, default is standard thermostat
      if( word == test ){
         if(value == "standard"){
            sld::internal::thermostat =
               sld::internal::standard_thermostat;
            return true;
         }
         if(value == "sled"){
            sld::internal::thermostat = sld::internal::sled_thermostat;
            return true;
         }
         err::zexit("spin-lattice:thermostat must be standard or sled");
      }

      test = "electron-temperature"; // initial electron temperature for SLED thermostat
      if( word == test ){
         double temperature = vin::str_to_double(value);
         vin::check_for_valid_value(temperature, word, line, prefix, unit,
                                    "none", 0.0, 1.0e6, "input",
                                    "0 - 1,000,000 K");
         sld::internal::initial_electron_temperature = temperature;
         sld::internal::initial_electron_temperature_set = true;
         return true;
      }

      test = "electron-spin-coupling"; // G_es [W m^-3 K^-1] for SLED thermostat
      if( word == test ){
         double coupling = vin::str_to_double(value);
         vin::check_for_valid_value(coupling, word, line, prefix, unit,
                                    "none", 0.0, 1.0e40, "input",
                                    "0 - 1E40 W m^-3 K^-1");
         sld::internal::electron_spin_coupling = coupling;
         sld::internal::electron_spin_coupling_set = true;
         return true;
      }

      test = "electron-spin-coupling-dynamic"; // derive G_es from the instantaneous spin state using Ma et al. Eq. (33)
      if( word == test ){
         sld::internal::electron_spin_coupling_dynamic =
            vin::check_for_valid_bool(value, word, line, prefix, "input");
         return true;
      }

      test = "electron-phonon-coupling"; // G_ep [W m^-3 K^-1] for SLED thermostat
      if( word == test ){
         double coupling = vin::str_to_double(value);
         vin::check_for_valid_value(coupling, word, line, prefix, unit,
                                    "none", 0.0, 1.0e40, "input",
                                    "0 - 1E40 W m^-3 K^-1");
         sld::internal::electron_phonon_coupling = coupling;
         return true;
      }

      test = "electron-heat-capacity-model"; // constant, linear or non-linear electron heat capacity model for SLED thermostat
      if( word == test ){
         if(value == "constant"){
            sld::internal::electron_heat_capacity_model =
               sld::internal::constant_electron_heat_capacity;
            return true;
         }
         if(value == "linear"){
            sld::internal::electron_heat_capacity_model =
               sld::internal::linear_electron_heat_capacity;
            return true;
         }
         if(value == "non-linear"){
            sld::internal::electron_heat_capacity_model =
               sld::internal::nonlinear_electron_heat_capacity;
            return true;
         }
         err::zexit("spin-lattice:electron-heat-capacity-model must be constant, linear or non-linear");
      }

      test = "electron-heat-capacity"; // C_e [J m^-3 K^-1] for SLED thermostat, only used if electron-heat-capacity-model = constant
      if( word == test ){
         double heat_capacity = vin::str_to_double(value);
         vin::check_for_valid_positive_value(
            heat_capacity, word, line, prefix, unit, "none", 1.0e-30,
            1.0e40, "input", "greater than zero in J m^-3 K^-1");
         sld::internal::electron_heat_capacity = heat_capacity;
         return true;
      }

      test = "electron-heat-capacity-coefficient"; // gamma [J m^-3 K^-2] for SLED thermostat, only used if electron-heat-capacity-model = linear
      if( word == test ){
         double coefficient = vin::str_to_double(value);
         vin::check_for_valid_positive_value(
            coefficient, word, line, prefix, unit, "none", 1.0e-30,
            1.0e40, "input", "greater than zero in J m^-3 K^-2");
         sld::internal::electron_heat_capacity_coefficient = coefficient;
         return true;
      }

      test = "potential-cutoff-range";
      if( word == test ){
         double r_c = vin::str_to_double(value);
         vin::check_for_valid_value(r_c, word, line, prefix, unit, "length", 2.0, 20.0,"input","2 - 20 A");
         sld::internal::r_cut_pot= r_c;
         return true;
      }

      test = "fields-cutoff-range";
      if( word == test ){
         double r_cf = vin::str_to_double(value);
         vin::check_for_valid_value(r_cf, word, line, prefix, unit, "length", 2, 20.0,"input","2 - 20 A");
         sld::internal::r_cut_fields= r_cf;
         if(!sld::internal::r_cut_exchange_set) sld::internal::r_cut_exchange = r_cf; // default to same cutoff for exchange if not set
         if(!sld::internal::r_cut_neel_l_set) sld::internal::r_cut_neel_l = r_cf;
         if(!sld::internal::r_cut_neel_q_set) sld::internal::r_cut_neel_q = r_cf;
         return true;
      }

      test = "exchange-cutoff-range";
      if( word == test ){
         double cutoff = vin::str_to_double(value);
         vin::check_for_valid_value(cutoff, word, line, prefix, unit, "length", 0.1, 20.0, "input", "0.1 - 20 A");
         sld::internal::r_cut_exchange = cutoff;
         sld::internal::r_cut_exchange_set = true;
         return true;
      }
      /* test = "fixed-lattice";
      if( word == test ){
      sld::internal::fixed-lattice=true;
      return true;
   }
   test = "fixed-spin";
   if( word == test ){
   sld::internal::fixed-spin=true;
   return true;
}*/

test = "initial-random-displacement";
if( word == test ){
   double dr_in = vin::str_to_double(value);
   vin::check_for_valid_value(dr_in, word, line, prefix, unit, "length", 0.001, 1.0,"input","0.001 - 1A");
   sld::internal::dr_init= dr_in;
   return true;
}

test = "initial-thermal-velocity";
if( word == test ){
   double temp = vin::str_to_double(value);
   vin::check_for_valid_value(temp, word, line, prefix, unit, "none", 0, 2000,"input","0 - 2000");
   sld::internal::th_velo= temp;
   return true;
}

//--------------------------------------------------------------------
// Keyword not found
//--------------------------------------------------------------------
return false;

}

//---------------------------------------------------------------------------
// Function to process material parameters
//---------------------------------------------------------------------------
bool match_material_parameter(std::string const word, std::string const value, std::string const unit, int const line, int const super_index, const int sub_index){

   // add prefix string
   std::string prefix="material:";

   // Check for material id > current array size and if so dynamically expand mp array
   if((unsigned int) super_index + 1 > internal::mp.size() && super_index + 1 < 101) internal::mp.resize(super_index + 1);
   std::string test = "mass";
   if( word == test ){
      double m = vin::str_to_double(value);
      vin::check_for_valid_value(m, word, line, prefix, unit, "mass", 1.0e-20, 1.0e20,"input","1E-20 - 1E20");
      sld::internal::mp[super_index].mass.set(m);
      return true;
   }

   test = "damping-constant-lattice";
   if( word == test ){
      double damp= vin::str_to_double(value);
      vin::check_for_valid_value(damp, word, line, prefix, unit, "none", 0, 1.0,"input","0- 1");
      sld::internal::mp[super_index].damp_lat.set(damp);
      return true;
   }

   test = "equilibration-damping-constant-lattice";
   if( word == test ){
      double damp= vin::str_to_double(value);
      vin::check_for_valid_value(damp, word, line, prefix, unit, "none", 0, 1.0,"input","0- 1");
      sld::internal::mp[super_index].eq_damp_lat.set(damp);
      return true;
   }

   test = "exchange-J0";
   if( word == test ){
      double j0 = vin::str_to_double(value);
      vin::check_for_valid_value(j0, word, line, prefix, unit, "energy", 0, 5,"input","0 - 5 eV");
      sld::internal::mp[super_index].J0.set(j0);
      return true;
   }

   test = "exchange-K0";
   if( word == test ){
      double k0 = vin::str_to_double(value);
      vin::check_for_valid_value(k0, word, line, prefix, unit, "energy", 0, 5,"input","+/- 5 eV");
      sld::internal::mp[super_index].K0.set(k0);
      return true;
   }

   test = "bethe-slater-alpha-j";
   if( word == test ){
      double parameter = vin::str_to_double(value);
      vin::check_for_valid_value(parameter, word, line, prefix, unit, "energy", 0, 5, "material", "+/- 5 eV");
      sld::internal::mp[super_index].bethe_slater_alpha_j.set(parameter);
      return true;
   }

   test = "bethe-slater-gamma-j";
   if( word == test ){
      double parameter = vin::str_to_double(value);
      vin::check_for_valid_value(parameter, word, line, prefix, unit, "none", 0, 100, "material", "+/- 100");
      sld::internal::mp[super_index].bethe_slater_gamma_j.set(parameter);
      return true;
   }

   test = "bethe-slater-delta-j";
   if( word == test ){
      double parameter = vin::str_to_double(value);
      vin::check_for_valid_positive_value(parameter, word, line, prefix, unit, "length", 0.01, 100, "material", "0.01 - 100 A");
      sld::internal::mp[super_index].bethe_slater_delta_j.set(parameter);
      return true;
   }

   test = "bethe-slater-alpha-k";
   if( word == test ){
      double parameter = vin::str_to_double(value);
      vin::check_for_valid_value(parameter, word, line, prefix, unit, "energy", 0, 5, "material", "+/- 5 eV");
      sld::internal::mp[super_index].bethe_slater_alpha_k.set(parameter);
      return true;
   }

   test = "bethe-slater-gamma-k";
   if( word == test ){
      double parameter = vin::str_to_double(value);
      vin::check_for_valid_value(parameter, word, line, prefix, unit, "none", 0, 100, "material", "+/- 100");
      sld::internal::mp[super_index].bethe_slater_gamma_k.set(parameter);
      return true;
   }

   test = "bethe-slater-delta-k";
   if( word == test ){
      double parameter = vin::str_to_double(value);
      vin::check_for_valid_positive_value(parameter, word, line, prefix, unit, "length", 0.01, 100, "material", "0.01 - 100 A");
      sld::internal::mp[super_index].bethe_slater_delta_k.set(parameter);
      return true;
   }

   test = "harmonic-potential-V0";
   if( word == test ){
      double v0 = vin::str_to_double(value);
      vin::check_for_valid_value(v0, word, line, prefix, unit, "energy", 1.0e-20, 1.0e20,"input","1E-20 - 1E20");
      sld::internal::mp[super_index].V0.set(v0);
      return true;
   }

   test = "coupling-C0";
   if( word == test ){
      double c0 = vin::str_to_double(value);
      vin::check_for_valid_value(c0, word, line, prefix, unit, "mass", 0, 1,"input","0 - 1");
      sld::internal::mp[super_index].C0.set(c0);
      return true;
   }

   test = "neel-C-l";
   if( word == test ){
      double coefficient = vin::str_to_double(value);
      vin::check_for_valid_value(coefficient, word, line, prefix, unit, "none", -10.0, 10.0, "material", "-10 - 10");
      sld::internal::mp[super_index].neel_C_l.set(coefficient);
      return true;
   }

   test = "neel-C-q";
   if( word == test ){
      double coefficient = vin::str_to_double(value);
      vin::check_for_valid_value(coefficient, word, line, prefix, unit, "none", -10.0, 10.0, "material", "-10 - 10");
      sld::internal::mp[super_index].neel_C_q.set(coefficient);
      return true;
   }

   test = "neel-alpha-l";
   if( word == test ){
      double parameter = vin::str_to_double(value);
      vin::check_for_valid_value(parameter, word, line, prefix, unit, "energy", -5.0, 5.0, "material", "+/- 5 eV");
      sld::internal::mp[super_index].neel_alpha_l.set(parameter);
      return true;
   }

   test = "neel-gamma-l";
   if( word == test ){
      double parameter = vin::str_to_double(value);
      vin::check_for_valid_value(parameter, word, line, prefix, unit, "none", -100.0, 100.0, "material", "+/- 100");
      sld::internal::mp[super_index].neel_gamma_l.set(parameter);
      return true;
   }

   test = "neel-delta-l";
   if( word == test ){
      double parameter = vin::str_to_double(value);
      vin::check_for_valid_positive_value(parameter, word, line, prefix, unit, "length", 0.01, 100.0, "material", "0.01 - 100 A");
      sld::internal::mp[super_index].neel_delta_l.set(parameter);
      return true;
   }

   test = "neel-alpha-q";
   if( word == test ){
      double parameter = vin::str_to_double(value);
      vin::check_for_valid_value(parameter, word, line, prefix, unit, "energy", -5.0, 5.0, "material", "+/- 5 eV");
      sld::internal::mp[super_index].neel_alpha_q.set(parameter);
      return true;
   }

   test = "neel-gamma-q";
   if( word == test ){
      double parameter = vin::str_to_double(value);
      vin::check_for_valid_value(parameter, word, line, prefix, unit, "none", -100.0, 100.0, "material", "+/- 100");
      sld::internal::mp[super_index].neel_gamma_q.set(parameter);
      return true;
   }

   test = "neel-delta-q";
   if( word == test ){
      double parameter = vin::str_to_double(value);
      vin::check_for_valid_positive_value(parameter, word, line, prefix, unit, "length", 0.01, 100.0, "material", "0.01 - 100 A");
      sld::internal::mp[super_index].neel_delta_q.set(parameter);
      return true;
   }


   //--------------------------------------------------------------------
   // Keyword not found
   //--------------------------------------------------------------------
   return false;

}



} // end of sld namespace
