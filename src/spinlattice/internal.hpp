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

#ifndef SLD_INTERNAL_H_
#define SLD_INTERNAL_H_
//
//---------------------------------------------------------------------
// This header file defines shared internal data structures and
// functions for the sld module. These functions and
// variables should not be accessed outside of this module.
//---------------------------------------------------------------------

// C++ standard library headers
#include <cmath>
#include <string>
#include <vector>

// Vampire headers
#include "sld.hpp"

// sld module headers
#include "internal.hpp"
#include "snap.hpp"


namespace sld{

   namespace internal{

      struct radial_result_t{
         double value;
         double derivative; // df/dr
      };

      struct neel_vector_t{
         double x;
         double y;
         double z;
      };

      struct neel_pair_result_t{
         double energy; // Phi_ij before the minus sign in H_Neel
         neel_vector_t spin_derivative_i; // dPhi_ij/ds_i
         neel_vector_t force_i; // dPhi_ij/dr_i in energy/Angstrom
         double spin_hessian_trace_i; // Tr[(I-s_i s_i^T) d^2 Phi_ij/ds_i^2]
      };

      // add two Cartesian vectors 
      inline neel_vector_t neel_add(const neel_vector_t& a, const neel_vector_t& b){
         neel_vector_t result = {a.x+b.x, a.y+b.y, a.z+b.z};
         return result;
      }

      // multiply a Cartesian vector by a scalar
      inline neel_vector_t neel_scale(const neel_vector_t& vector, const double scale){
         neel_vector_t result = {scale*vector.x, scale*vector.y, scale*vector.z};
         return result;
      }

      // calculate the Cartesian dot product of two vectors
      inline double neel_dot(const neel_vector_t& a, const neel_vector_t& b){
         return a.x*b.x+a.y*b.y+a.z*b.z;
      }

      // calculate f(r) and df/dr for the Bethe-Slater radial function
      inline radial_result_t bethe_slater(const double alpha,
                                          const double gamma,
                                          const double delta,
                                          const double r){
         const double inverse_delta_squared = 1.0/(delta*delta);
         const double x = r*r*inverse_delta_squared;
         const double exponential = std::exp(-x);
         radial_result_t result;
         result.value = 4.0*alpha*x*(1.0-gamma*x)*exponential;
         result.derivative = 8.0*alpha*r*inverse_delta_squared*exponential*
                             (1.0-x-gamma*x*(2.0-x));
         return result;
      }

      // calculate f(r) and df/dr for the inverse fourth radial function
      inline radial_result_t inverse_fourth(const double coefficient,
                                            const double r){
         const double inverse_r = 1.0/r;
         const double inverse_r2 = inverse_r*inverse_r;
         radial_result_t result;
         result.value = coefficient*inverse_r2*inverse_r2;
         result.derivative = -4.0*result.value*inverse_r;
         return result;
      }

      // apply a hard radial cutoff to a value and its derivative
      inline radial_result_t apply_hard_cutoff(const radial_result_t& radial,
                                                const double r,
                                                const double cutoff){
         if(r < cutoff) return radial;
         radial_result_t zero = {0.0, 0.0};
         return zero;
      }

      // apply a smooth cutoff between the switch and cutoff distances
      inline radial_result_t apply_smooth_cutoff(const radial_result_t& radial,
                                                 const double r,
                                                 const double switch_distance,
                                                 const double cutoff){
         if(r <= switch_distance) return radial;
         if(r >= cutoff){
            radial_result_t zero = {0.0, 0.0};
            return zero;
         }

         const double inverse_width = 1.0/(cutoff-switch_distance);
         const double x = (r-switch_distance)*inverse_width;
         const double x2 = x*x;
         const double x3 = x2*x;
         const double x4 = x3*x;
         const double x5 = x4*x;
         const double switching = 1.0-10.0*x3+15.0*x4-6.0*x5;
         const double switching_derivative = (-30.0*x2+60.0*x3-30.0*x4)*inverse_width;

         radial_result_t result;
         // d[f(r)S(r)]/dr=f'(r)S(r)+f(r)S'(r)
         result.value = radial.value*switching;
         result.derivative = radial.derivative*switching+
                             radial.value*switching_derivative;
         return result;
      }

      // calculate the energy, spin derivative, force and curvature of one full Neel pair
      inline neel_pair_result_t evaluate_neel_pair(const neel_vector_t& displacement,
                                                   const neel_vector_t& spin_i,
                                                   const neel_vector_t& spin_j,
                                                   const radial_result_t& l_curve,
                                                   const radial_result_t& q_curve){

         neel_pair_result_t result;
         result.energy = 0.0;
         result.spin_derivative_i = neel_vector_t{0.0, 0.0, 0.0};
         result.force_i = neel_vector_t{0.0, 0.0, 0.0};
         result.spin_hessian_trace_i = 0.0;

         const double distance_squared = neel_dot(displacement, displacement);
         const double distance = std::sqrt(distance_squared);
         const double inverse_distance = 1.0/distance;
         const neel_vector_t bond = neel_scale(displacement, inverse_distance);

         const double x = neel_dot(bond, spin_i);
         const double y = neel_dot(bond, spin_j);
         const double z = neel_dot(spin_i, spin_j);
         const double x2 = x*x;
         const double y2 = y*y;
         const double x3 = x2*x;
         const double y3 = y2*y;

         // A=xy-z/3, B=(x^2-z/3)(y^2-z/3), C=xy^3+yx^3
         const double one_third = 1.0/3.0;
         const double A = x*y-one_third*z;
         const double u = x2-one_third*z;
         const double v = y2-one_third*z;
         const double B = u*v;
         const double C = x*y3+y*x3;

         // l1=l+12q/35, q1=9q/5 and q2=-2q/5
         const double l1 = l_curve.value+(12.0/35.0)*q_curve.value;
         const double q1 = (9.0/5.0)*q_curve.value;
         const double q2 = -(2.0/5.0)*q_curve.value;
         const double l1_derivative = l_curve.derivative+(12.0/35.0)*q_curve.derivative;
         const double q1_derivative = (9.0/5.0)*q_curve.derivative;
         const double q2_derivative = -(2.0/5.0)*q_curve.derivative;

         // Phi_ij=l*A+(12/35)q*A+(9/5)q*B-(2/5)q*C=l1(r)A+q1(r)B+q2(r)C
         result.energy = l1*A+q1*B+q2*C;

         // dPhi/ds_i=l1*dA/ds_i+q1*dB/ds_i+q2*dC/ds_i
         const neel_vector_t dA_dsi = neel_add(neel_scale(bond, y),
                                               neel_scale(spin_j, -one_third));
         const neel_vector_t du_dsi = neel_add(neel_scale(bond, 2.0*x),
                                               neel_scale(spin_j, -one_third));
         const neel_vector_t dv_dsi = neel_scale(spin_j, -one_third);
         const neel_vector_t dB_dsi = neel_add(neel_scale(du_dsi, v),
                                               neel_scale(dv_dsi, u));
         const neel_vector_t dC_dsi =
            neel_scale(bond, y3+3.0*y*x2);
         result.spin_derivative_i = neel_add(neel_scale(dA_dsi, l1),
                                    neel_add(neel_scale(dB_dsi, q1),
                                             neel_scale(dC_dsi, q2)));

         // tau=Tr[(I-s_i s_i^T)d^2 Phi/ds_i^2].  The l1*A term is linear in s_i and has zero curvature  
         // for the nonlinear terms tau_B=2y^2-2x^2y^2+2x^2z-4xy/3-2z/3+2(1-z^2)/9 and tau_C=6xy(1-x^2)
         const double tau_B = 2.0*y2-2.0*x2*y2+2.0*x2*z
                            -(4.0/3.0)*x*y-(2.0/3.0)*z
                            +(2.0/9.0)*(1.0-z*z);
         const double tau_C = 6.0*x*y*(1.0-x2);
         result.spin_hessian_trace_i = q1*tau_B+q2*tau_C;

         // grad(x)=(s_i-x*e)/r and grad(y)=(s_j-y*e)/r
         const neel_vector_t grad_x =
            neel_scale(neel_add(spin_i, neel_scale(bond, -x)),
                       inverse_distance);
         const neel_vector_t grad_y =
            neel_scale(neel_add(spin_j, neel_scale(bond, -y)),
                       inverse_distance);
         const neel_vector_t grad_A =
            neel_add(neel_scale(grad_x, y), neel_scale(grad_y, x));
         const neel_vector_t grad_B =
            neel_add(neel_scale(grad_x, 2.0*x*v),
                     neel_scale(grad_y, 2.0*y*u));
         const neel_vector_t grad_C =
            neel_add(neel_scale(grad_x, y3+3.0*y*x2),
                     neel_scale(grad_y, 3.0*x*y2+x3));

         // F_i=grad(Phi)=Phi'(r)e+l1*grad(A)+q1*grad(B)+q2*grad(C)
         const double radial_derivative =
            l1_derivative*A+q1_derivative*B+q2_derivative*C;
         result.force_i =  neel_add(neel_scale(bond, radial_derivative),
                           neel_add(neel_scale(grad_A, l1),
                                    neel_add(neel_scale(grad_B, q1),
                                             neel_scale(grad_C, q2))));
         return result;
      }

      class set_double_t{

      private:
         double value; // value
         bool setf; // flag specifiying variable has been set

      public:
         // class functions
         // constructor
         set_double_t() : value(0.0), setf(false) { }

         // setting function
         void set(double in_value){
            value = in_value;
            setf = true;
         };

         // get value function
         double get(){ return value; };
         // check if variable is set
         bool is_set(){ return setf; };

      };

      //-------------------------------------------------------------------------
      // Internal data type definitions
      //-------------------------------------------------------------------------

      //-----------------------------------------------------------------------------
      // internal materials class for storing material parameters
      //-----------------------------------------------------------------------------
      class mp_t{

          private:

          public:

             //------------------------------
             // material parameter variables
             //------------------------------
             set_double_t mass;
             set_double_t V0;
             set_double_t J0;
             set_double_t C0;
             set_double_t neel_C_l;
             set_double_t neel_C_q;
             set_double_t damp_lat;
             set_double_t eq_damp_lat;

             set_double_t J0_ms;
             set_double_t K0;
             set_double_t K0_ms;
             set_double_t C0_ms;
             set_double_t J0_prime;
             set_double_t K0_prime;
             set_double_t bethe_slater_alpha_j;
             set_double_t bethe_slater_gamma_j;
             set_double_t bethe_slater_delta_j;
             set_double_t bethe_slater_alpha_k;
             set_double_t bethe_slater_gamma_k;
             set_double_t bethe_slater_delta_k;
             set_double_t neel_alpha_l;
             set_double_t neel_gamma_l;
             set_double_t neel_delta_l;
             set_double_t neel_alpha_q;
             set_double_t neel_gamma_q;
             set_double_t neel_delta_q;
             set_double_t F_th_sigma;
             set_double_t F_th_sigma_eq;







             // constructor
             mp_t (const unsigned int max_materials = 100){
                mass.set(5.7915e-3);
                V0.set(0.15);
                J0.set(0.904);
                J0_prime.set(3*0.904/7.8);
                J0_ms.set(0.904/2.04028e-23);
                C0.set(0.5);
                C0_ms.set(0.5/2.04028e-23);
                neel_C_l.set(1.0e-2); 
                neel_C_q.set(-1.0e-4);
                F_th_sigma.set(1.0);
                F_th_sigma_eq.set(1.0);
                damp_lat.set(0.06);
                eq_damp_lat.set(0.6);



             }; // end of constructor

             

       }; // end of internal::mp class


       // NEW: Parameters for the linear phonon pump
      extern bool linear_pump_enabled;
      extern double phonon_frequency;
      extern double phonon_force_amplitude[3]; // For fx, fy, fz
      extern double phonon_wavevector[3];      // For kx, ky, kz
      extern double phonon_pulse_start_time;
      extern double phonon_pulse_end_time;
      extern std::vector<double> coupling_field_x; 
      extern std::vector<double> coupling_field_y; 
      extern std::vector<double> coupling_field_z;
      // sum of Tr[(I-ss^T)dB/ds] in tesla over the supported nonlinear hamiltonian terms used for spin temperature calculations
      extern std::vector<double> spin_hessian_trace;

      // NEW: Intermediate storage for wave parameters from the input file
      extern double phonon_wave_lambda[3];   // For λ_x, λ_y, λ_z
      extern double phonon_wave_direction[3]; // For direction_x, direction_y, direction_z
 


      //-------------------------------------------------------------------------
      // Internal shared variables
      //-------------------------------------------------------------------------

      enum lattice_potential_t{
         no_lattice_potential,
         harmonic_lattice_potential,
         morse_lattice_potential,
         snap_lattice_potential,
         snap_zbl_lattice_potential
      };

      enum exchange_function_t{
         cubic_exchange_function,
         bethe_slater_exchange_function
      };

      enum spin_hamiltonian_t{
         bilinear_spin_hamiltonian,
         biquadratic_spin_hamiltonian
      };

      enum neel_radial_function_t{
         inverse_fourth_neel_radial_function,
         bethe_slater_neel_radial_function
      };

      enum neel_cutoff_function_t{
         hard_neel_cutoff_function,
         smooth_neel_cutoff_function
      };

      enum thermostat_t{
         standard_thermostat,
         sled_thermostat
      };

      enum electron_heat_capacity_t{
         constant_electron_heat_capacity,
         linear_electron_heat_capacity,
         nonlinear_electron_heat_capacity
      };

      extern bool enabled; // bool to enable module
      extern std::vector<sld::internal::mp_t> mp; // array of material properties

      extern double r_cut_pot; // mechanical potential cutoff
      extern double r_cut_fields; // maximum exchange/coupling neighbour cutoff
      extern double r_cut_exchange;
      extern double r_cut_neel_l;
      extern double r_cut_neel_q;
      extern double r_switch_neel_l;
      extern double r_switch_neel_q;
      extern bool r_cut_exchange_set;
      extern bool r_cut_neel_l_set;
      extern bool r_cut_neel_q_set;
      extern bool r_switch_neel_l_set;
      extern bool r_switch_neel_q_set;
      extern exchange_function_t exchange_function;
      extern spin_hamiltonian_t spin_hamiltonian;
      extern neel_radial_function_t neel_radial_function;
      extern neel_cutoff_function_t neel_cutoff_function;
      extern bool exchange_offset;
      extern thermostat_t thermostat;
      extern electron_heat_capacity_t electron_heat_capacity_model;
      extern double initial_electron_temperature;
      extern bool initial_electron_temperature_set;
      extern double electron_temperature;
      extern double electron_spin_coupling;
      extern double electron_phonon_coupling;
      extern bool electron_spin_coupling_dynamic;
      extern bool electron_spin_coupling_set;
      extern double electron_heat_capacity;
      extern double electron_heat_capacity_coefficient;
      extern bool sled_production_initialized;

      extern double dr_init;
      extern double th_velo;

       //for the morse potential
       extern double morse_beta;
       extern double morse_factor;
       extern double alpha_m;
       extern double r0_m;
       extern double morse_D;

      extern lattice_potential_t lattice_potential;
      extern bool pseudodipolar;
      extern bool full_neel;
      extern bool harmonic_debug_enabled;
      extern int harmonic_debug_force_calls;
      extern int harmonic_debug_max_force_calls;
      extern double zbl_inner_cutoff;
      extern double zbl_outer_cutoff;
      extern double zbl_atomic_number;
      extern snap_potential_t snap_potential;

      bool lattice_potential_is_mlip();

      void prepare_sled_thermostat();
      void update_sled_thermostat();
      void initialise_sled_couplings();
      double get_electron_heat_capacity(const double temperature);

      void print_force_debug_summary(const std::string& label,
                                     const int call,
                                     const int start_index,
                                     const int end_index,
                                     const std::vector<int>& neighbour_counts,
                                     const std::vector<double>& initial_forces_x,
                                     const std::vector<double>& initial_forces_y,
                                     const std::vector<double>& initial_forces_z,
                                     const std::vector<double>& forces_array_x,
                                     const std::vector<double>& forces_array_y,
                                     const std::vector<double>& forces_array_z);


      //extern std::vector<int> sld_neighbour_list_start_index;
      //extern std::vector<int> sld_neighbour_list_end_index;
      //extern std::vector<int> sld_neighbour_list_array;

      extern std::vector<double> x0_coord_array;
      extern std::vector<double> y0_coord_array;
      extern std::vector<double> z0_coord_array;


      extern std::vector <double> x_coord_storage_array;
      extern std::vector <double> y_coord_storage_array;
      extern std::vector <double> z_coord_storage_array;



      extern std::vector<double> forces_array_x;
      extern std::vector<double> forces_array_y;
      extern std::vector<double> forces_array_z;

      extern std::vector<double> fields_array_x;
      extern std::vector<double> fields_array_y;
      extern std::vector<double> fields_array_z;

      extern std::vector<double> velo_array_x;
      extern std::vector<double> velo_array_y;
      extern std::vector<double> velo_array_z;
      extern std::vector<double> potential_eng;
      extern std::vector<double> exch_eng;
      extern std::vector<double> coupl_eng;
      extern std::vector<double> sumJ;
      extern std::vector<double> sumC;

      extern std::vector<int> test_atom_list; //Core atoms of each octant



      void initialise_positions(std::vector<double>& x0_coord_array, // coord vectors for atoms
                  std::vector<double>& y0_coord_array,
                  std::vector<double>& z0_coord_array,
                  std::vector<double>& x_coord_array, // coord vectors for atoms
                  std::vector<double>& y_coord_array,
                  std::vector<double>& z_coord_array,
                  const double dr);

      void thermal_velocity(std::vector<double>& x_velo_array, // coord vectors for atoms
                 std::vector<double>& y_velo_array,
                 std::vector<double>& z_velo_array);

//function to resize vectors and initialise rest of parameters
      void initialise_sld_parameters();

//functions to compute potentials
      void compute_forces_harmonic(const int start_index,
            const int end_index, // last +1 atom to be calculated
            const std::vector<int>& neighbour_list_start_index,
            const std::vector<int>& neighbour_list_end_index,
            const std::vector<int>& type_array, // type for atom
            const std::vector<int>& neighbour_list_array, // list of interactions between atom
            const std::vector<double>& x0_coord_array, // coord vectors for atoms
            const std::vector<double>& y0_coord_array,
            const std::vector<double>& z0_coord_array,
            const std::vector<double>& x_coord_array, // coord vectors for atoms
            const std::vector<double>& y_coord_array,
            const std::vector<double>& z_coord_array,
            std::vector<double>& forces_array_x, //  vectors for forces
            std::vector<double>& forces_array_y,
            std::vector<double>& forces_array_z,
            std::vector<double>& potential_eng);

//functions to compute fields
      void compute_exchange(const int start_index,
            const int end_index, // last +1 atom to be calculated
            const std::vector<int>& neighbour_list_start_index,
            const std::vector<int>& neighbour_list_end_index,
            const std::vector<int>& type_array, // type for atom
            const std::vector<int>& neighbour_list_array, // list of interactions between atom
            const std::vector<double>& x_coord_array, // coord vectors for atoms
            const std::vector<double>& y_coord_array,
            const std::vector<double>& z_coord_array,
            const std::vector<double>& x_spin_array, // coord vectors for atoms
            const std::vector<double>& y_spin_array,
            const std::vector<double>& z_spin_array,
            std::vector<double>& forces_array_x, //  vectors for forces
            std::vector<double>& forces_array_y,
            std::vector<double>& forces_array_z,
            std::vector<double>& fields_array_x, //  vectors for forces
            std::vector<double>& fields_array_y,
            std::vector<double>& fields_array_z);


      void compute_forces_morse(const int start_index,
            const int end_index, // last +1 atom to be calculated
            const std::vector<int>& neighbour_list_start_index,
            const std::vector<int>& neighbour_list_end_index,
            const std::vector<int>& type_array, // type for atom
            const std::vector<int>& neighbour_list_array, // list of interactions between atom
            const std::vector<double>& x_coord_array, // coord vectors for atoms
            const std::vector<double>& y_coord_array,
            const std::vector<double>& z_coord_array,
            std::vector<double>& forces_array_x, //  vectors for forces
            std::vector<double>& forces_array_y,
            std::vector<double>& forces_array_z,
            std::vector<double>& potential_eng);

      void compute_thz(const int start_index,
            const int end_index, 
            const std::vector<double>& x_coord_array, // current coord vectors for atoms
            const std::vector<double>& y_coord_array,
            const std::vector<double>& z_coord_array,
            std::vector<double>& forces_array_x, //  vectors for forces
            std::vector<double>& forces_array_y,
            std::vector<double>& forces_array_z);

      void compute_forces_snap(const int start_index,
            const int end_index, // last +1 atom to be calculated
            const std::vector<int>& neighbour_list_start_index,
            const std::vector<int>& neighbour_list_end_index,
            const std::vector<int>& type_array, // type for atom
            const std::vector<int>& neighbour_list_array, // list of interactions between atom
            const std::vector<double>& x_coord_array, // coord vectors for atoms
            const std::vector<double>& y_coord_array,
            const std::vector<double>& z_coord_array,
            std::vector<double>& forces_array_x, //  vectors for forces
            std::vector<double>& forces_array_y,
            std::vector<double>& forces_array_z,
            std::vector<double>& potential_eng);

      void compute_forces_snap_zbl(const int start_index,
            const int end_index, // last +1 atom to be calculated
            const std::vector<int>& neighbour_list_start_index,
            const std::vector<int>& neighbour_list_end_index,
            const std::vector<int>& type_array, // type for atom
            const std::vector<int>& neighbour_list_array, // list of interactions between atom
            const std::vector<double>& x_coord_array, // coord vectors for atoms
            const std::vector<double>& y_coord_array,
            const std::vector<double>& z_coord_array,
            std::vector<double>& forces_array_x, //  vectors for forces
            std::vector<double>& forces_array_y,
            std::vector<double>& forces_array_z,
            std::vector<double>& potential_eng);
//
      void compute_sld_coupling(const int start_index,
            const int end_index, // last +1 atom to be calculated
            const std::vector<int>& neighbour_list_start_index,
            const std::vector<int>& neighbour_list_end_index,
            const std::vector<int>& type_array, // type for atom
            const std::vector<int>& neighbour_list_array, // list of interactions between atom
            const std::vector<double>& x_coord_array, // coord vectors for atoms
            const std::vector<double>& y_coord_array,
            const std::vector<double>& z_coord_array,
            const std::vector<double>& x_spin_array, // spin  vectors for atoms
            const std::vector<double>& y_spin_array,
            const std::vector<double>& z_spin_array,
            std::vector<double>& forces_array_x, //  vectors for forces
            std::vector<double>& forces_array_y,
            std::vector<double>& forces_array_z,
            std::vector<double>& fields_array_x, //  vectors for forces
            std::vector<double>& fields_array_y,
            std::vector<double>& fields_array_z);

      // Calculate the full neel coupling fields, forces, energies and spin curvature
      void compute_sld_coupling_neel(const int start_index,
            const int end_index, // last +1 atom to be calculated
            const std::vector<int>& neighbour_list_start_index,
            const std::vector<int>& neighbour_list_end_index,
            const std::vector<int>& type_array, // type for atom
            const std::vector<int>& neighbour_list_array, // list of interactions between atom
            const std::vector<double>& x_coord_array, // coord vectors for atoms
            const std::vector<double>& y_coord_array,
            const std::vector<double>& z_coord_array,
            const std::vector<double>& x_spin_array, // spin  vectors for atoms
            const std::vector<double>& y_spin_array,
            const std::vector<double>& z_spin_array,
            std::vector<double>& forces_array_x, //  vectors for forces
            std::vector<double>& forces_array_y,
            std::vector<double>& forces_array_z,
            std::vector<double>& fields_array_x, //  vectors for forces
            std::vector<double>& fields_array_y,
            std::vector<double>& fields_array_z);

//

      void cayley_update(const int start_index,
                  const int end_index,
                  double dt,
                  std::vector<double>& x_spin_array, // coord vectors for atoms
                  std::vector<double>& y_spin_array,
                  std::vector<double>& z_spin_array,
                  std::vector<double>& fields_array_x, //  vectors for fields
                  std::vector<double>& fields_array_y,
                  std::vector<double>& fields_array_z);
   // add spin noise

      void add_spin_noise(const int start_index,
                  const int end_index,
                  double dt,
                  const std::vector<int>& type_array, // type for atom
                  const std::vector<double>& x_spin_array, // coord vectors for atoms
                  const std::vector<double>& y_spin_array,
                  const std::vector<double>& z_spin_array,
                  std::vector<double>& fields_array_x, //  vectors for fields
                  std::vector<double>& fields_array_y,
                  std::vector<double>& fields_array_z,
                  std::vector<double>& Hx_th, //  vectors for fields
                  std::vector<double>& Hy_th,
                  std::vector<double>& Hz_th);

    //MPI variables
    extern std::vector<std::vector<int> > c_octants; //Core atoms of each octant
    extern std::vector<std::vector<int> > b_octants; //Boundary atoms of each octant

    extern std::vector <int> all_atoms_octant_start_index;
    extern std::vector <int> all_atoms_octant_end_index;
    extern std::vector <int> all_atoms_octant;


      //-------------------------------------------------------------------------
      // Internal function declarations
      //-------------------------------------------------------------------------

   } // end of internal namespace

} // end of sld namespace


#endif //SLD_INTERNAL_H_
