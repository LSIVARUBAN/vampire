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

   // restore the SLED thermostat state after continuing from a checkpoint
   void restore_thermostat_checkpoint(){
      if(internal::thermostat != internal::sled_thermostat && internal::thermostat != internal::sled_energy_thermostat) return;
      // sim::temperature is T_e so the checkpointed temperature is the recorded T_e
      internal::electron_temperature = sim::temperature;
      internal::sled_production_initialized = sim::time > sim::equilibration_time; 
   }

namespace internal{

namespace{

   double sled_volume_m3 = 0.0;
   double electron_number_density = 0.0;

   // return the fixed number density used by the SLED heat capacity
   double get_electron_number_density(){

      if(electron_number_density != 0.0) 
         return electron_number_density;

      double number_of_atoms = 0.0;
      #ifdef MPICF
         number_of_atoms = static_cast<double>(vmpi::num_core_atoms+vmpi::num_bdry_atoms);
         MPI_Allreduce(MPI_IN_PLACE, &number_of_atoms, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #else
         number_of_atoms = static_cast<double>(atoms::num_atoms);
      #endif

      electron_number_density = number_of_atoms/sled_volume_m3;
      return electron_number_density;
   }

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
            // Non-linear Fe heat capacity from Ma et al.,
            // Phys. Rev. B 85, 184301 (2012):
            // https://doi.org/10.1103/PhysRevB.85.184301
            return 3.0 * constants::kB * get_electron_number_density() * std::tanh(2.0e-4 * temperature);
         }
      }

      return 0.0;
   }

namespace{

   // return u_e(T)=integral_0^T C_e(T')dT' in J m^-3
   double get_electron_energy_density(const double temperature){

      switch(electron_heat_capacity_model){
         case constant_electron_heat_capacity:
            return electron_heat_capacity * temperature;

         case linear_electron_heat_capacity:
            return 0.5 * electron_heat_capacity_coefficient * temperature * temperature;

         case nonlinear_electron_heat_capacity:{
            const double heat_capacity_limit = 3.0 * constants::kB * get_electron_number_density();
            const double coefficient = 2.0e-4;
            const double scaled_temperature = coefficient * temperature;
            double log_cosh = 0.0;
            if(scaled_temperature < 20.0){
               log_cosh = std::log(std::cosh(scaled_temperature));
            }
            else{
               log_cosh = scaled_temperature - std::log(2.0) + std::log1p(std::exp(-2.0*scaled_temperature));
            }
            return heat_capacity_limit/coefficient * log_cosh;
         }
      }

      return 0.0;
   }

   // solve the electron energy equation for T_e at the end of the step
   double update_electron_temperature_from_energy(const double temperature, const double spin_temperature, const double lattice_temperature){

      const double total_coupling = electron_spin_coupling + electron_phonon_coupling;
      if(total_coupling == 0.0) 
         return temperature;

      const double target_temperature = (electron_spin_coupling*spin_temperature + electron_phonon_coupling*lattice_temperature) / total_coupling;
      if(target_temperature == temperature) 
         return temperature;

      // u_e(T_e^{n+1})+dt(G_es+G_ep)T_e^{n+1}
      // =u_e(T_e^n)+dt(G_es T_s+G_ep T_l)
      const double right_hand_side = get_electron_energy_density(temperature) + mp::dt_SI * total_coupling * target_temperature;

      double lower_temperature = std::min(temperature, target_temperature);
      double upper_temperature = std::max(temperature, target_temperature);
      double new_temperature = 0.5*(lower_temperature+upper_temperature);

      for(int iteration = 0; iteration < 32; ++iteration){
         const double residual = get_electron_energy_density(new_temperature) + mp::dt_SI * total_coupling * new_temperature - right_hand_side;

         if(residual == 0.0) 
            return new_temperature;

         if(residual > 0.0) 
            upper_temperature = new_temperature;
         else lower_temperature = new_temperature;

         if(upper_temperature-lower_temperature <= 1.0e-12*std::max(1.0, upper_temperature)){
            return 0.5*(lower_temperature+upper_temperature);
         }

         const double derivative = get_electron_heat_capacity(new_temperature) + mp::dt_SI * total_coupling;
         const double newton_temperature = new_temperature-residual/derivative;

         if(newton_temperature > lower_temperature && newton_temperature < upper_temperature){
            if(std::abs(newton_temperature-new_temperature) <= 1.0e-12*std::max(1.0, new_temperature)){
               return newton_temperature;
            }
            new_temperature = newton_temperature;
         }
         else{
            new_temperature = 0.5*(lower_temperature+upper_temperature);
         }
      }

      return new_temperature;
   }

} // end of anonymous namespace

   // prepare the SLED thermostat for the current time step, setting the electron temperature to the equilibration temperature during equilibration and to the initial electron temperature at the start of production
   void prepare_sled_thermostat(){

      if(thermostat != sled_thermostat && thermostat != sled_energy_thermostat) return;

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

      if(thermostat != sled_thermostat && thermostat != sled_energy_thermostat) return;

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
         electron_spin_coupling = get_dynamic_electron_spin_coupling(0, end_index);
      }

      if(thermostat == sled_energy_thermostat){
         electron_temperature = update_electron_temperature_from_energy(electron_temperature, sld::spin_temperature, sld::lattice_temperature);
      }
      else{
         // retain the original explicit Euler SLED temperature update
         const double heat_capacity = get_electron_heat_capacity(electron_temperature);
         const double energy_transfer = electron_spin_coupling * (electron_temperature-sld::spin_temperature) + electron_phonon_coupling * (electron_temperature-sld::lattice_temperature);
         electron_temperature -= mp::dt_SI * energy_transfer/heat_capacity;
         electron_temperature = std::max(0.0, electron_temperature);
      }
      sim::temperature = electron_temperature;
   }

} // namespace internal
} // namespace sld
