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
#include <vector>

// Vampire headers
#include "sld.hpp"
#include "material.hpp"

// sld module headers
#include "internal.hpp"

namespace sld{

   //------------------------------------------------------------------------------
   // Externally visible variables
   //------------------------------------------------------------------------------
   bool enabled = false;
   double var_test=0;
   double spin_temperature;
   double lattice_temperature;
   double J_eff;
   double C_eff;
   bool suzuki_trotter_parallel_initialized = false;

   namespace internal{

      //------------------------------------------------------------------------
      // Shared variables inside sld module
      //------------------------------------------------------------------------

      //Spin-lattice coupling variables
      bool linear_pump_enabled;
      double phonon_frequency;
      double phonon_force_amplitude[3];
      double phonon_wavevector[3];
      double phonon_pulse_start_time;
      double phonon_pulse_end_time;
      double phonon_wave_lambda[3];   
      double phonon_wave_direction[3]; 
      std::vector<double> coupling_field_x; 
      std::vector<double> coupling_field_y; 
      std::vector<double> coupling_field_z;
      std::vector<double> spin_hessian_trace;
      bool spin_temperature_correction = true;


      bool enabled; // bool to enable module

      std::vector<internal::mp_t> mp; // array of material properties

      double r_cut_pot; // mechanical potential cutoff
      double r_cut_fields = 3.75; // default cutoff for exchange and Neel fields
      double r_cut_exchange = 3.75;
      double r_cut_neel_l = 3.75;
      double r_cut_neel_q = 3.75;
      double r_switch_neel_l = 0.0;
      double r_switch_neel_q = 0.0;
      bool r_cut_exchange_set = false;
      bool r_cut_neel_l_set = false;
      bool r_cut_neel_q_set = false;
      bool r_switch_neel_l_set = false;
      bool r_switch_neel_q_set = false;
      exchange_function_t exchange_function = cubic_exchange_function; // default exchange function
      spin_hamiltonian_t spin_hamiltonian = bilinear_spin_hamiltonian; // default spin hamiltonian
      neel_radial_function_t neel_radial_function = inverse_fourth_neel_radial_function;
      neel_cutoff_function_t neel_cutoff_function = hard_neel_cutoff_function;
      bool exchange_offset = false; // by default don't offset the hamiltonian
      thermostat_t thermostat = standard_thermostat;
      electron_heat_capacity_t electron_heat_capacity_model = linear_electron_heat_capacity;
      double initial_electron_temperature = 300.0; // T_e0 [K]
      bool initial_electron_temperature_set = false;
      double electron_temperature = 300.0; // T_e [K]
      double electron_spin_coupling = 0.0; // G_es [W m^-3 K^-1]
      double electron_phonon_coupling = 0.0; // G_ep [W m^-3 K^-1]
      bool electron_spin_coupling_dynamic = false; // derive G_es from Ma et al. Eq. (33)
      bool electron_spin_coupling_set = false;
      double electron_heat_capacity = 1.0; // [J m^-3 K^-1]
      double electron_heat_capacity_coefficient = 225.0; // [J m^-3 K^-2]
      bool sled_production_initialized = false;

      double dr_init; // initial conditions
      double th_velo = 0.0;

      double morse_beta;
      double morse_factor;
      double alpha_m;
      double r0_m;
      double morse_D;

      lattice_potential_t lattice_potential = no_lattice_potential;
      bool pseudodipolar;
      bool full_neel;
      bool harmonic_debug_enabled = false; // ouput force info
      int harmonic_debug_force_calls = 0; // number of calls to debug forces
      int harmonic_debug_max_force_calls = 5;
      double zbl_inner_cutoff = 4.0; // default cutoff for ZBL (bcc Fe)
      double zbl_outer_cutoff = 4.8;
      double zbl_atomic_number = 26.0;
      snap_potential_t snap_potential; // snap potential class

      bool lattice_potential_is_mlip(){

         switch(lattice_potential){
            case snap_lattice_potential:
            case snap_zbl_lattice_potential:
               return true;
            default:
               return false;
         }

      }


      //initial sld neighbor list
      //std::vector<int> sld_neighbour_list_start_index;
      //std::vector<int> sld_neighbour_list_end_index;
      //std::vector<int> sld_neighbour_list_array;

      std::vector<double> x0_coord_array;
      std::vector<double> y0_coord_array;
      std::vector<double> z0_coord_array;


      std::vector <double> x_coord_storage_array;
      std::vector <double> y_coord_storage_array;
      std::vector <double> z_coord_storage_array;

      std::vector<double> forces_array_x;
      std::vector<double> forces_array_y;
      std::vector<double> forces_array_z;

      std::vector<double> fields_array_x;
      std::vector<double> fields_array_y;
      std::vector<double> fields_array_z;

      std::vector<double> velo_array_x;
      std::vector<double> velo_array_y;
      std::vector<double> velo_array_z;

      std::vector<double> potential_eng;
      std::vector<double> sumJ;
      std::vector<double> sumC;
      std::vector<double> exch_eng;
      std::vector<double> coupl_eng;

      std::vector<int> test_atom_list; //Core atoms of each octant

      //MPI variables
      std::vector<std::vector<int> > c_octants; //Core atoms of each octant
      std::vector<std::vector<int> > b_octants; //Boundary atoms of each octant
      std::vector <int> all_atoms_octant_start_index;
      std::vector <int> all_atoms_octant_end_index;
      std::vector <int> all_atoms_octant;


   } // end of internal namespace

   // nonlinear field curvature array
   const std::vector<double>& spin_temperature_hessian_trace(){
      return internal::spin_hessian_trace;
   }

   // return whether nonlinear spin temperature corrections are enabled
   bool spin_temperature_correction_enabled(){
      return internal::spin_temperature_correction;
   }

} // end of sld namespace
