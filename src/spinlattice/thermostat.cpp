//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   Dev: Lennon Sivaruban, 2026.
//
//   Email: lennon.sivaruban@postgrad.manchester.ac.uk
//
//------------------------------------------------------------------------------

// C++ standard library headers
#include <algorithm>
#include <cmath>

// Vampire headers
#include "atoms.hpp"
#include "constants.hpp"
#include "create.hpp"
#include "errors.hpp"
#include "material.hpp"
#include "sim.hpp"
#include "sld.hpp"
#include "vmpi.hpp"

// sld module headers
#include "internal.hpp"

namespace sld{

   // return the current electron temperature
   double electron_temperature(){
      return internal::electron_temperature;
   }

namespace internal{

namespace{

   double sled_volume_m3 = 0.0;

   // Ma et al., Phys. Rev. B 85, 184301 (2012), Eq. (33):
   // G_es = 2 k_B/(hbar V) sum_i gamma_es,i <s_i . H_i>.
   // In VAMPIRE units this becomes
   // G_es = 2 k_B gamma/V sum_i [alpha_i/(1+alpha_i^2)] s_i.H_i.
   double get_dynamic_electron_spin_coupling(const int start_index,
                                              const int end_index){

      double weighted_spin_field = 0.0;
      for(int atom = start_index; atom < end_index; ++atom){
         const int material = atoms::type_array[atom];
         const double alpha = ::mp::material[material].alpha;
         const double spin_dot_field =
              atoms::x_spin_array[atom] * fields_array_x[atom]
            + atoms::y_spin_array[atom] * fields_array_y[atom]
            + atoms::z_spin_array[atom] * fields_array_z[atom];
         weighted_spin_field +=
            alpha / (1.0 + alpha * alpha) * spin_dot_field;
      }

      #ifdef MPICF
         MPI_Allreduce(MPI_IN_PLACE, &weighted_spin_field, 1,
                       MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif

      const double coupling =
         2.0 * constants::kB * ::mp::gamma_SI * weighted_spin_field /
         sled_volume_m3;

      return std::max(0.0, coupling); // keep above zero
   }

} // end of anonymous namespace

   // Validate the SLED coupling inputs
   void initialise_sled_couplings(){

      if(electron_spin_coupling_dynamic && electron_spin_coupling_set){
         err::zexit("Specify either spin-lattice:electron-spin-coupling or spin-lattice:electron-spin-coupling-dynamic");
      }
      // coupling coefficients and heat capacities are volume densities... use the same simulation box volume as the SLED heat capacity model
      sled_volume_m3 =
         cs::system_dimensions[0] * cs::system_dimensions[1] *
         cs::system_dimensions[2] * 1.0e-30;
   }

   // return the current electron heat capacity
   double get_electron_heat_capacity(const double temperature){

      switch(electron_heat_capacity_model){
         case constant_electron_heat_capacity:
            return electron_heat_capacity;

         case linear_electron_heat_capacity:
            return electron_heat_capacity_coefficient * temperature;

         case nonlinear_electron_heat_capacity:{
            // NV simulations so just calculate the density once
            static double number_density = 0.0;
            if(number_density == 0.0){
               double number_of_atoms = 0.0;
               #ifdef MPICF
                  number_of_atoms = static_cast<double>(vmpi::num_core_atoms+vmpi::num_bdry_atoms);
                  MPI_Allreduce(MPI_IN_PLACE, &number_of_atoms, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
               #else
                  number_of_atoms = static_cast<double>(atoms::num_atoms);
               #endif

               const double volume_m3 =
                  cs::system_dimensions[0] * cs::system_dimensions[1] * cs::system_dimensions[2] * 1.0e-30;
               number_density = number_of_atoms/volume_m3;
            }

            // Non-linear Fe heat capacity from Ma et al.,
            // Phys. Rev. B 85, 184301 (2012):
            // https://doi.org/10.1103/PhysRevB.85.184301
            return 3.0 * constants::kB * number_density * std::tanh(2.0e-4 * temperature);
         }
      }

      return 0.0;
   }

   // prepare the SLED thermostat for the current time step, setting the electron temperature to the equilibration temperature during equilibration and to the initial electron temperature at the start of production
   void prepare_sled_thermostat(){

      if(thermostat != sled_thermostat) return;

      // set the electron temperature to the equilibration temperature during equilibration
      if(sim::time < sim::equilibration_time){
         electron_temperature = sim::Teq;
         sled_production_initialized = false;
      }
      // set the electron temperature to the initial electron temperature at the start of production
      else if(!sled_production_initialized){
         electron_temperature = initial_electron_temperature;
         sled_production_initialized = true;
      }

      // set the Langevin bath temperature to the current electron temperature
      sim::temperature = electron_temperature;
   }

   // update the SLED thermostat, calculating the new electron temperature based on the current spin and lattice temperatures
   void update_sled_thermostat(){

      if(thermostat != sled_thermostat) return;

      const int end_index =
         #ifdef MPICF
            vmpi::num_core_atoms+vmpi::num_bdry_atoms;
         #else
            atoms::num_atoms;
         #endif

      // reset the fields and forces arrays to zero before calculating the spin temperature
      std::fill(fields_array_x.begin(), fields_array_x.end(), 0.0);
      std::fill(fields_array_y.begin(), fields_array_y.end(), 0.0);
      std::fill(fields_array_z.begin(), fields_array_z.end(), 0.0);
      std::fill(forces_array_x.begin(), forces_array_x.end(), 0.0);
      std::fill(forces_array_y.begin(), forces_array_y.end(), 0.0);
      std::fill(forces_array_z.begin(), forces_array_z.end(), 0.0);

      // rebuild the fields and forces arrays and calculate the current spin and lattice temperatures
      sld::compute_fields(0,
                          end_index,
                          atoms::neighbour_list_start_index,
                          atoms::neighbour_list_end_index,
                          atoms::type_array,
                          atoms::neighbour_list_array,
                          atoms::x_coord_array,
                          atoms::y_coord_array,
                          atoms::z_coord_array,
                          atoms::x_spin_array,
                          atoms::y_spin_array,
                          atoms::z_spin_array,
                          forces_array_x,
                          forces_array_y,
                          forces_array_z,
                          fields_array_x,
                          fields_array_y,
                          fields_array_z);

      sld::spin_temperature = sld::compute_spin_temperature(0,
                                       end_index,
                                       atoms::type_array,
                                       atoms::x_spin_array,
                                       atoms::y_spin_array,
                                       atoms::z_spin_array,
                                       fields_array_x,
                                       fields_array_y,
                                       fields_array_z,
                                       mp::mu_s_array);

      sld::lattice_temperature = sld::compute_lattice_temperature(0,
                                          end_index,
                                          atoms::type_array,
                                          atoms::x_velo_array,
                                          atoms::y_velo_array,
                                          atoms::z_velo_array);

      // equilibration phase uses the shared bath. energy exchange starts in production 
      if(sim::time < sim::equilibration_time){
         electron_temperature = sim::Teq;
         sim::temperature = electron_temperature;
         return;
      }

      if(electron_spin_coupling_dynamic){
         electron_spin_coupling =
            get_dynamic_electron_spin_coupling(0, end_index); 
      }

      // calculate the heat capacity
      const double heat_capacity = get_electron_heat_capacity(electron_temperature);
      // calculate the energy transfer between the electron, spin and lattice systems based on the current temperatures and coupling constants
      const double energy_transfer = electron_spin_coupling * (electron_temperature - sld::spin_temperature) 
                                 + electron_phonon_coupling * (electron_temperature - sld::lattice_temperature);

      // update the electron temperature
      electron_temperature -= mp::dt_SI * energy_transfer / heat_capacity; 
      electron_temperature = std::max(0.0, electron_temperature);
      sim::temperature = electron_temperature;
   }

} // namespace internal
} // namespace sld
