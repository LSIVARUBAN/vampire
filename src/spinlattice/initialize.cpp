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
#include <algorithm>
#include <cmath>
#include <sstream>

// Vampire headers
#include "sld.hpp"
#include "atoms.hpp"
#include "iostream"
#include "neighbours.hpp"
#include "material.hpp"
#include "constants.hpp"
#include "errors.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "vmpi.hpp"

// sld module headers
#include "internal.hpp"

namespace sld{

   //----------------------------------------------------------------------------
   // Function to initialize sld module
   //----------------------------------------------------------------------------
   void initialize(){

      // check for sld module being enabled
      if(!sld::enabled) return;

      std::cout<<"Input parameters for Spin-lattice dynamics simulations:"<<std::endl;
      std::cout<<"*******************************************************"<<std::endl;

      // MLIP potential files define the cutoff, rather than the VAMPIRE input file potential-cutoff-range
      if(sld::internal::lattice_potential == sld::internal::snap_lattice_potential ||
         sld::internal::lattice_potential == sld::internal::snap_zbl_lattice_potential){
         sld::internal::snap_potential.initialise(mp::num_materials);
         sld::internal::r_cut_pot = sld::internal::snap_potential.cutoff(); // assign snap cutoff to global var
      }
      if(sld::internal::lattice_potential == sld::internal::snap_zbl_lattice_potential){
         if(sld::internal::zbl_inner_cutoff >= sld::internal::zbl_outer_cutoff){
            err::zexit("spin-lattice:zbl-inner-cutoff must be smaller than spin-lattice:zbl-outer-cutoff"); // make sure zbl cutoffs are appropriate
         }
         sld::internal::r_cut_pot = std::max(sld::internal::r_cut_pot,
                                             sld::internal::zbl_outer_cutoff);
      }
      if((sld::internal::lattice_potential == sld::internal::snap_lattice_potential ||
          sld::internal::lattice_potential == sld::internal::snap_zbl_lattice_potential) &&
         sld::internal::r_cut_pot > sld::internal::r_cut_fields){
         // SNAP filters the exchange derived neighbour list, so its complete mechanical cutoff must fit inside the range available to that list
         err::zexit("SNAP potential cutoff exceeds spin-lattice:fields-cutoff-range");
      }

      // check that exchange and neel cutoffs are smaller than the fields cutoff
      if(sld::internal::r_cut_exchange > sld::internal::r_cut_fields){
         err::zexit("spin-lattice:exchange-cutoff-range exceeds spin-lattice:fields-cutoff-range");
      }
      if(sld::internal::full_neel && (sld::internal::r_cut_neel_l > sld::internal::r_cut_fields || sld::internal::r_cut_neel_q > sld::internal::r_cut_fields)){
         err::zexit("Full-Neel cutoff exceeds spin-lattice:fields-cutoff-range");
      }

      const bool use_smooth_neel_cutoff = sld::internal::neel_cutoff_function == sld::internal::smooth_neel_cutoff_function;
      if(sld::internal::full_neel && use_smooth_neel_cutoff){
         if(!sld::internal::r_switch_neel_l_set ||
            !sld::internal::r_switch_neel_q_set){
            err::zexit("Smooth full Neel cutoff needs neel-l and neel-q ranges to be set");
         }
         if(sld::internal::r_switch_neel_l >= sld::internal::r_cut_neel_l ||
            sld::internal::r_switch_neel_q >= sld::internal::r_cut_neel_q){
            err::zexit("The full Neel range must be smaller than its cutoff range");
         }
      }

      std::cout<<"Potential Cutoff: "<<sld::internal::r_cut_pot<<std::endl;
      std::cout<<"Maximum fields cutoff: "<<sld::internal::r_cut_fields<<std::endl;
      std::cout<<"Exchange cutoff: "<<sld::internal::r_cut_exchange<<std::endl;
      std::cout<<"Mass: "<<sld::internal::mp[0].mass.get()<<std::endl;
      std::cout<<"Lattice damping: "<<sld::internal::mp[0].damp_lat.get()<<std::endl;
      std::cout<<"Coupling C0 "<<sld::internal::mp[0].C0.get()<<std::endl;

      const bool use_bethe_slater =
         sld::internal::exchange_function ==
         sld::internal::bethe_slater_exchange_function;
      const bool use_biquadratic =
         sld::internal::spin_hamiltonian ==
         sld::internal::biquadratic_spin_hamiltonian;

      // validate each curve as a complete alpha/gamma/delta triple
      if(use_bethe_slater){
         for(int mat = 0; mat < mp::num_materials; ++mat){
            const bool j_parameters_set =
               sld::internal::mp[mat].bethe_slater_alpha_j.is_set() &&
               sld::internal::mp[mat].bethe_slater_gamma_j.is_set() &&
               sld::internal::mp[mat].bethe_slater_delta_j.is_set();
            if(!j_parameters_set){
               err::zexit("Bethe-Slater exchange requires alpha-j, gamma-j and delta-j in every material");
            }

            const bool k_parameters_set =
               sld::internal::mp[mat].bethe_slater_alpha_k.is_set() &&
               sld::internal::mp[mat].bethe_slater_gamma_k.is_set() &&
               sld::internal::mp[mat].bethe_slater_delta_k.is_set();
            if(use_biquadratic && !k_parameters_set){
               err::zexit("Bethe-Slater biquadratic exchange requires alpha-k, gamma-k and delta-k in every material");
            }
         }
      }
      else if(use_biquadratic){
         for(int mat = 0; mat < mp::num_materials; ++mat){
            if(!sld::internal::mp[mat].K0.is_set()){
               err::zexit("Cubic biquadratic exchange requires exchange-K0 in every material");
            }
         }
      }

      const bool use_bethe_slater_neel = sld::internal::neel_radial_function == sld::internal::bethe_slater_neel_radial_function;
      if(sld::internal::full_neel && use_bethe_slater_neel){
         for(int mat = 0; mat < mp::num_materials; ++mat){
            const bool l_parameters_set =
               sld::internal::mp[mat].neel_alpha_l.is_set() &&
               sld::internal::mp[mat].neel_gamma_l.is_set() &&
               sld::internal::mp[mat].neel_delta_l.is_set();
            const bool q_parameters_set =
               sld::internal::mp[mat].neel_alpha_q.is_set() &&
               sld::internal::mp[mat].neel_gamma_q.is_set() &&
               sld::internal::mp[mat].neel_delta_q.is_set();
            if(!l_parameters_set || !q_parameters_set){
               err::zexit("Bethe-Slater full Neel coupling needs complete alpha/gamma/delta sets for l and q in every material");
            }
         }
      }
      std::cout << "Exchange function: " << (use_bethe_slater ? "Bethe-Slater" : "cubic") << std::endl;
      std::cout << "Spin Hamiltonian: " << (use_biquadratic ? "biquadratic" : "bilinear") << std::endl;
      std::cout << "Exchange Hamiltonian offset: " << (sld::internal::exchange_offset ? "enabled" : "disabled") << std::endl;

      const bool use_sled_thermostat = sld::internal::thermostat == sld::internal::sled_thermostat || sld::internal::thermostat == sld::internal::sled_energy_thermostat;
      if(use_sled_thermostat){
         if(!sld::internal::initial_electron_temperature_set){
            // if the initial electron temperature is not set, use the equilibration temperature as the initial value
            sld::internal::initial_electron_temperature = sim::temperature;
         }
         sld::internal::electron_temperature = sim::Teq; 
         sld::internal::sled_production_initialized = false;
         sld::internal::initialise_sled_couplings();

         const double initial_heat_capacity = sld::internal::get_electron_heat_capacity(sld::internal::initial_electron_temperature);
         if(initial_heat_capacity <= 0.0){
            err::zexit("The SLED electron heat capacity must be greater than zero");
         }

         std::cout << "Thermostat: " << (sld::internal::thermostat == sld::internal::sled_energy_thermostat ? "SLED energy" : "SLED") << std::endl;
         std::cout << "Initial electron temperature: " << sld::internal::initial_electron_temperature << " K" << std::endl;
         if(sld::internal::electron_spin_coupling_dynamic){
            std::cout << "Electron-spin coupling (G_es): dynamic" << std::endl;
         }
         else{
            std::cout << "Electron-spin coupling (G_es): " << sld::internal::electron_spin_coupling << " W m^-3 K^-1" << std::endl;
         }
         std::cout << "Electron-phonon coupling (G_ep): " << sld::internal::electron_phonon_coupling << " W m^-3 K^-1" << std::endl;
      }
      else{
         std::cout << "Thermostat: standard" << std::endl;
      }

      switch(sld::internal::lattice_potential){
         case sld::internal::harmonic_lattice_potential:
            std::cout<<"Harmonic potential is used of potential well depth V0="<<sld::internal::mp[0].V0.get()<<std::endl;
            if(sld::internal::harmonic_debug_enabled){
               std::cout<<"Harmonic debug output enabled for the first " <<sld::internal::harmonic_debug_max_force_calls <<" force calls"<<std::endl;
            }
            break;
         case sld::internal::morse_lattice_potential:
            std::cout<<"Morse potential is used"<<std::endl;
            break;
         case sld::internal::snap_lattice_potential:
            std::cout<<"SNAP potential is used with "<<sld::internal::snap_potential.number_of_elements()<<" element(s) and "<<sld::internal::snap_potential.number_of_coefficients()<<" descriptor coefficients"<<std::endl;
            break;
         case sld::internal::snap_zbl_lattice_potential:
            std::cout<<"SNAP+ZBL potential is used with "<<sld::internal::snap_potential.number_of_elements()<<" element(s) and "<<sld::internal::snap_potential.number_of_coefficients()<<" descriptor coefficients"<<std::endl;
            std::cout<<"ZBL atomic number: "<<sld::internal::zbl_atomic_number<<std::endl;
            std::cout<<"ZBL inner cutoff: "<<sld::internal::zbl_inner_cutoff<<" Angstrom"<<std::endl;
            std::cout<<"ZBL outer cutoff: "<<sld::internal::zbl_outer_cutoff<<" Angstrom"<<std::endl;
            break;
         default:
            std::cout<<"No lattice potential selected"<<std::endl;
            break;
      }

      if(sld::internal::pseudodipolar)std::cout<<"Pseudodipolar coupling is used of strength C0="<<sld::internal::mp[0].C0.get()<<std::endl;
      if(sld::internal::full_neel){
         std::cout<<"Full Neel coupling is used with " << (use_bethe_slater_neel ? "Bethe-Slater" : "inverse-fourth") << " radial functions and " << (use_smooth_neel_cutoff ? "smooth" : "hard") <<" cutoffs"<<std::endl;
         std::cout<<"Full Neel l cutoff: "<<sld::internal::r_cut_neel_l << " Angstrom"<<std::endl;
         std::cout<<"Full Neel q cutoff: "<<sld::internal::r_cut_neel_q <<" Angstrom"<<std::endl;
         if(!use_bethe_slater_neel){
            std::cout<<"Full Neel C_l: "<<sld::internal::mp[0].neel_C_l.get()<<std::endl;
            std::cout<<"Full Neel C_q: "<<sld::internal::mp[0].neel_C_q.get()<<std::endl;
         }
      }

      std::cout<<"*******************************************************"<<std::endl;


      sld::internal::coupling_field_x.resize(atoms::num_atoms, 0.0);
      sld::internal::coupling_field_y.resize(atoms::num_atoms, 0.0);
      sld::internal::coupling_field_z.resize(atoms::num_atoms, 0.0);
      sld::internal::spin_hessian_trace.resize(atoms::num_atoms, 0.0); // for spin temperature calculation

      // Author: Muhammad Hamza Asim
      // Calculate phonon wavevector components (k) from wavelength and direction
      if (sld::internal::linear_pump_enabled) {

          // Define constants for the calculation
          const double epsilon = 1.0e-12; // A small number to prevent division by zero
          const double two_pi = 6.28318530718;

          
          // The formula is k = (2*pi / lambda) * direction_vector.
          if (std::abs(sld::internal::phonon_wave_lambda[0]) > epsilon) {
              sld::internal::phonon_wavevector[0] = two_pi / sld::internal::phonon_wave_lambda[0] * sld::internal::phonon_wave_direction[0];
          } else {
              sld::internal::phonon_wavevector[0] = 0.0;
          }

          if (std::abs(sld::internal::phonon_wave_lambda[1]) > epsilon) {
              sld::internal::phonon_wavevector[1] = two_pi / sld::internal::phonon_wave_lambda[1] * sld::internal::phonon_wave_direction[1];
          } else {
              sld::internal::phonon_wavevector[1] = 0.0;
          }

          if (std::abs(sld::internal::phonon_wave_lambda[2]) > epsilon) {
              sld::internal::phonon_wavevector[2] = two_pi / sld::internal::phonon_wave_lambda[2] * sld::internal::phonon_wave_direction[2];
          } else {
              sld::internal::phonon_wavevector[2] = 0.0;
          }



      }


    //initialise exchange, coupling parameters
     sld::internal::initialise_sld_parameters();

     //for the morse potential:

     sld::internal::alpha_m= 1.3885;
     sld::internal::r0_m=2.845;
     sld::internal::morse_D=0.4174;;
     sld::internal::morse_beta=exp( sld::internal::alpha_m * sld::internal::r0_m);
     sld::internal::morse_factor=-2.0 *sld::internal::morse_D * sld::internal::alpha_m;

     //initialise for Parallel simulations
     //Initialize parallel mc variables
     suzuki_trotter_parallel_initialized = false;
     internal::c_octants.resize(8);
     internal::b_octants.resize(8);

     internal::all_atoms_octant_start_index.reserve(8);
     internal::all_atoms_octant_end_index.reserve(8);
     internal::all_atoms_octant.reserve(atoms::num_atoms);


     sld::internal::initialise_positions(sld::internal::x0_coord_array, // coord vectors for atoms
                  sld::internal::y0_coord_array,
                  sld::internal::z0_coord_array,
                  atoms::x_coord_array, // coord vectors for atoms
                  atoms::y_coord_array,
                  atoms::z_coord_array,
                  sld::internal::dr_init);

     if(sld::internal::th_velo > 0.0){ // initialise thermal velocities
        sld::internal::thermal_velocity(atoms::x_velo_array,
                                        atoms::y_velo_array,
                                        atoms::z_velo_array);
     }

    // sld::tests();



      return;


   }


   namespace internal{

   void initialise_positions(std::vector<double>& x0_coord_array, // coord vectors for atoms
               std::vector<double>& y0_coord_array,
               std::vector<double>& z0_coord_array,
               std::vector<double>& x_coord_array, // coord vectors for atoms
               std::vector<double>& y_coord_array,
               std::vector<double>& z_coord_array,
               const double dr){


   x0_coord_array.resize(atoms::num_atoms,0);
   y0_coord_array.resize(atoms::num_atoms,0);
   z0_coord_array.resize(atoms::num_atoms,0);



   x_coord_storage_array.resize(atoms::num_atoms,0);
   y_coord_storage_array.resize(atoms::num_atoms,0);
   z_coord_storage_array.resize(atoms::num_atoms,0);

   for( int i = 0; i < atoms::num_atoms; i++){

      x0_coord_array[i]=x_coord_array[i];
      y0_coord_array[i]=y_coord_array[i];
      z0_coord_array[i]=z_coord_array[i];

     x_coord_storage_array[i]=x_coord_array[i];
     y_coord_storage_array[i]=y_coord_array[i];
     z_coord_storage_array[i]=z_coord_array[i];



   }
   //std::cout<<"random positions  "<<dr<<std::endl;
   /*
    for( int i = 0; i < atoms::num_atoms; i++)
    {

             x_coord_array[i] += dr* (2.0*rand()/double(RAND_MAX) -1.0);
             y_coord_array[i] += dr* (2.0*rand()/double(RAND_MAX) -1.0);
             z_coord_array[i] += dr* (2.0*rand()/double(RAND_MAX) -1.0);


    }*/


   return;


   }//end of initialise initialise_positions


   // generate mass dependent thermal velocities with zero com motion for faster equilibration
   void thermal_velocity(std::vector<double>& x_velo_array, // velocities vectors
                         std::vector<double>& y_velo_array,
                         std::vector<double>& z_velo_array){
      int number_of_local_atoms = atoms::num_atoms;
      #ifdef MPICF
         number_of_local_atoms = vmpi::num_core_atoms+vmpi::num_bdry_atoms;
      #endif

      double total_mass = 0.0;
      double momentum_x = 0.0;
      double momentum_y = 0.0;
      double momentum_z = 0.0;
      double number_of_atoms = 0.0;

      // <v_alpha^2>=k_B T/m using velocities in Angstrom/ps
      for(int atom=0; atom<number_of_local_atoms; ++atom){
         const unsigned int material = atoms::type_array[atom];
         const double mass = sld::internal::mp[material].mass.get();
         const double sigma = std::sqrt(constants::kB_eV*sld::internal::th_velo/mass);
         x_velo_array[atom] = sigma*mtrandom::gaussian();
         y_velo_array[atom] = sigma*mtrandom::gaussian();
         z_velo_array[atom] = sigma*mtrandom::gaussian();

         total_mass += mass;
         momentum_x += mass*x_velo_array[atom];
         momentum_y += mass*y_velo_array[atom];
         momentum_z += mass*z_velo_array[atom];
         number_of_atoms += 1.0;
      }

      #ifdef MPICF
         double sums[5] = {
            total_mass, momentum_x, momentum_y, momentum_z, number_of_atoms
         };
         MPI_Allreduce(MPI_IN_PLACE, sums, 5, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         total_mass = sums[0];
         momentum_x = sums[1];
         momentum_y = sums[2];
         momentum_z = sums[3];
         number_of_atoms = sums[4];
      #endif

      if(total_mass <= 0.0 || number_of_atoms <= 1.0) return;

      const double centre_of_mass_velocity_x = momentum_x/total_mass;
      const double centre_of_mass_velocity_y = momentum_y/total_mass;
      const double centre_of_mass_velocity_z = momentum_z/total_mass;
      double thermal_mass_velocity_squared = 0.0;
      for(int atom=0; atom<number_of_local_atoms; ++atom){
         x_velo_array[atom] -= centre_of_mass_velocity_x;
         y_velo_array[atom] -= centre_of_mass_velocity_y;
         z_velo_array[atom] -= centre_of_mass_velocity_z;

         const unsigned int material = atoms::type_array[atom];
         const double mass = sld::internal::mp[material].mass.get();
         thermal_mass_velocity_squared += mass*(
            x_velo_array[atom]*x_velo_array[atom]+
            y_velo_array[atom]*y_velo_array[atom]+
            z_velo_array[atom]*z_velo_array[atom]);
      }

      #ifdef MPICF
         MPI_Allreduce(MPI_IN_PLACE, &thermal_mass_velocity_squared, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif

      // removing com translation leaves 3n-3 kinetic DOFs
      const double degrees_of_freedom = 3.0*number_of_atoms-3.0;
      const double scale = std::sqrt(
         degrees_of_freedom*constants::kB_eV*sld::internal::th_velo/
         thermal_mass_velocity_squared);
      for(int atom=0; atom<number_of_local_atoms; ++atom){
         x_velo_array[atom] *= scale;
         y_velo_array[atom] *= scale;
         z_velo_array[atom] *= scale;
      }
      return;

   }

   void initialise_sld_parameters(){

      sld::internal::forces_array_x.resize(atoms::num_atoms,0);
      sld::internal::forces_array_y.resize(atoms::num_atoms,0);
      sld::internal::forces_array_z.resize(atoms::num_atoms,0);

      sld::internal::fields_array_x.resize(atoms::num_atoms,0);
      sld::internal::fields_array_y.resize(atoms::num_atoms,0);
      sld::internal::fields_array_z.resize(atoms::num_atoms,0);

      sld::internal::velo_array_x.resize(atoms::num_atoms,0);
      sld::internal::velo_array_y.resize(atoms::num_atoms,0);
      sld::internal::velo_array_z.resize(atoms::num_atoms,0);

      // allocate atom indexed buffers once instead of on every integration step. Integrators overwrite each entry with new random nos before using them
      sld::internal::spin_noise_array_x.resize(atoms::num_atoms,0);
      sld::internal::spin_noise_array_y.resize(atoms::num_atoms,0);
      sld::internal::spin_noise_array_z.resize(atoms::num_atoms,0);
      sld::internal::lattice_noise_array_x.resize(atoms::num_atoms,0);
      sld::internal::lattice_noise_array_y.resize(atoms::num_atoms,0);
      sld::internal::lattice_noise_array_z.resize(atoms::num_atoms,0);

      // allocate one set of integration coefficients per material (removes sqrts, unit conversions and parameter selection from the atom loops)
      sld::internal::spin_damping_array.resize(mp::num_materials,0);
      sld::internal::spin_noise_scale_array.resize(mp::num_materials,0);
      sld::internal::lattice_damping_factor_array.resize(mp::num_materials,0);
      sld::internal::lattice_noise_scale_array.resize(mp::num_materials,0);
      sld::internal::lattice_dt2_over_mass_array.resize(mp::num_materials,0);

      sld::internal::potential_eng.resize(atoms::num_atoms,0);
      sld::internal::sumJ.resize(atoms::num_atoms,0);
      sld::internal::sumC.resize(atoms::num_atoms,0);
      sld::internal::exch_eng.resize(atoms::num_atoms,0);
      sld::internal::coupl_eng.resize(atoms::num_atoms,0);





      //sqrt( 2.0 * eta * consts::kB * T / (mass *dt));
      for(int mat=0;mat<mp::num_materials; mat++){

         double C=sld::internal::mp[mat].C0.get();
         sld::internal::mp[mat].C0.set(sld::internal::mp[mat].J0.get()*C);
         sld::internal::mp[mat].J0_ms.set(sld::internal::mp[mat].J0.get()/mp::material[mat].mu_s_SI);
         if(sld::internal::mp[mat].K0.is_set()){
            sld::internal::mp[mat].K0_ms.set(
               sld::internal::mp[mat].K0.get()/mp::material[mat].mu_s_SI);
            sld::internal::mp[mat].K0_prime.set(
               3.0*sld::internal::mp[mat].K0.get()/sld::internal::r_cut_exchange);
         }
         sld::internal::mp[mat].C0_ms.set(sld::internal::mp[mat].C0.get()/mp::material[mat].mu_s_SI);
         sld::internal::mp[mat].J0_prime.set(3.0*sld::internal::mp[mat].J0.get()/sld::internal::r_cut_exchange);
         sld::internal::mp[mat].F_th_sigma.set(sqrt(2.0*sld::internal::mp[mat].damp_lat.get()*constants::kB_eV / (sld::internal::mp[mat].mass.get()*mp::dt_SI*1e12)));
         sld::internal::mp[mat].F_th_sigma_eq.set(sqrt(2.0*sld::internal::mp[mat].eq_damp_lat.get()*constants::kB_eV / (sld::internal::mp[mat].mass.get()*mp::dt_SI*1e12)));
         //now change to ev the following
         sld::internal::mp[mat].V0.set(sld::internal::mp[mat].V0.get()*6.242e18);

         //initialise mass for the rest of the simulations
         mp::material[mat].mass=sld::internal::mp[mat].mass.get();
       }


       //initialise mass array


       for(int atom=0;atom<atoms::num_atoms;atom++){

   		int mat=atoms::type_array[atom];
        atoms::mass_spin_array[atom]=mp::material[mat].mass;
        }



   return;
   }

   // record repeated integrator coefficients for each material to avoid repeated caluclations in the integrator loops
   void prepare_integrator_coefficients(){

      const bool equilibrating = sim::time < sim::equilibration_time; // so we know which set of params to use
      const double sqrt_temperature = std::sqrt(sim::temperature); // both spin and lattice noise scale with sqrt(T)
      const double half_dt_ps = 0.5*mp::dt_SI*1.0e12; // lattice velocities are in ps and the integrator advances in half time steps

      // compute each value once per material since all atoms of the material share these params
      for(int mat = 0; mat < mp::num_materials; ++mat){
         spin_damping_array[mat] = equilibrating ? ::mp::material[mat].alpha_eq : ::mp::material[mat].alpha; // gilbert damping for current stage
         spin_noise_scale_array[mat] = sqrt_temperature*(equilibrating ? ::mp::material[mat].H_th_sigma_eq : ::mp::material[mat].H_th_sigma); // this scales the spin field to its amplitude for the current stage
         const double damping = equilibrating ? sld::internal::mp[mat].eq_damp_lat.get() : sld::internal::mp[mat].damp_lat.get(); // select the lattice damping for the current stage
         lattice_damping_factor_array[mat] = 1.0 - half_dt_ps*damping; // 1-gamma*dt/2
         lattice_noise_scale_array[mat] = sqrt_temperature*(equilibrating ? sld::internal::mp[mat].F_th_sigma_eq.get() : sld::internal::mp[mat].F_th_sigma.get()); // this scales the lattice force
         lattice_dt2_over_mass_array[mat] = half_dt_ps/sld::internal::mp[mat].mass.get(); // dt_2_m=dt/(2m) converts force into the velocity change for one half timestep
      }
   }

} //end of internal
} // end of sld namespace
