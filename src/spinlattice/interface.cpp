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

// Vampire headers
#include "sld.hpp"
#include "errors.hpp"
#include "vio.hpp"
#include "iostream"

// sld module headers
#include "internal.hpp"
#include "sld.hpp"


namespace sld{

   //---------------------------------------------------------------------------
   // Function to process input file parameters for sld module
   //---------------------------------------------------------------------------
   bool match_input_parameter(std::string const key, std::string const word, std::string const value, std::string const unit, int const line){

      // Check for valid key, if no match return false
      std::string prefix="spin-lattice";
      if(key!=prefix) return false;
      //------------------------------------------------------------------------
      // If spin-transport parameter is requested, then enable module
      //------------------------------------------------------------------------
      internal::enabled=true;
      //----------------------------------
      // Now test for all valid options
      //----------------------------------



      std::string test = "potential";
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
          return true;
         }
         test="full-neel";
         if( value == test ){
          sld::internal::full_neel=true;
          return true;
         }
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
     

      //--------------------------------------------------------------------
      // Keyword not found
      //--------------------------------------------------------------------
      return false;

   }



} // end of sld namespace
