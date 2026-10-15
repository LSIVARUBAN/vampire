//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Mara Strungaru 2022. All rights reserved.
//
//   Email: mara.strungaru@york.ac.uk
//
//   implementation based on the paper Phys. Rev. B 103, 024429, (2021) M.Strungaru, M.O.A. Ellis et al
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>

// Vampire headers
#include "anisotropy.hpp"
#include "atoms.hpp"
#include "create.hpp"
#include "material.hpp"
#include "sld.hpp"
#include "sim.hpp"

// sld module headers
#include "internal.hpp"

namespace sld{

   // calculate all spin fields and their mechanical forces
   void compute_fields(const int start_index, // first atom for exchange interactions to be calculated
               const int end_index,
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
               std::vector<double>& fields_array_x, //  vectors for fields
               std::vector<double>& fields_array_y,
               std::vector<double>& fields_array_z){

      // reset chi_i=Tr[(I-s_i s_i^T)dB_i/ds_i] before accumulating nonlinear terms
      std::fill(sld::internal::spin_hessian_trace.begin()+start_index,
                sld::internal::spin_hessian_trace.begin()+end_index,
                0.0);

      internal::compute_exchange(start_index, end_index,
                                 neighbour_list_start_index, neighbour_list_end_index,
                                 type_array, neighbour_list_array,
                                 x_coord_array, y_coord_array, z_coord_array,
                                 x_spin_array, y_spin_array, z_spin_array,
                                 forces_array_x, forces_array_y, forces_array_z,
                                 fields_array_x, fields_array_y, fields_array_z);

      if(sld::internal::pseudodipolar){ internal::compute_sld_coupling(start_index, end_index,
                                        neighbour_list_start_index, neighbour_list_end_index,
                                        type_array, neighbour_list_array,
                                        x_coord_array, y_coord_array, z_coord_array,
                                        x_spin_array, y_spin_array, z_spin_array,
                                        forces_array_x, forces_array_y, forces_array_z,
                                        fields_array_x, fields_array_y, fields_array_z);
      }

      if(sld::internal::full_neel){ internal::compute_sld_coupling_neel(start_index, end_index,
                                             neighbour_list_start_index, neighbour_list_end_index,
                                             type_array, neighbour_list_array,
                                             x_coord_array, y_coord_array, z_coord_array,
                                             x_spin_array, y_spin_array, z_spin_array,
                                             forces_array_x, forces_array_y, forces_array_z,
                                             fields_array_x, fields_array_y, fields_array_z);
      }

      internal::add_applied_field(start_index, end_index, fields_array_x, fields_array_y, fields_array_z);

      // add the anisotropy fields and their nonlinear curvature
      anisotropy::fields(atoms::x_spin_array, atoms::y_spin_array,
                         atoms::z_spin_array, atoms::type_array,
                         fields_array_x, fields_array_y, fields_array_z,
                         start_index, end_index, sim::temperature);
      anisotropy::spin_temperature_curvature(
         x_spin_array, y_spin_array, z_spin_array, type_array,
         sld::internal::spin_hessian_trace, start_index, end_index);

      return;
   }

namespace internal{

   // add the applied field after equilibration, the step from sim::time = equilibration_time is the first production step
   void add_applied_field(const int start_index,
                          const int end_index,
                          std::vector<double>& fields_array_x,
                          std::vector<double>& fields_array_y,
                          std::vector<double>& fields_array_z){
      if(sim::time >= sim::equilibration_time || sim::equilibration_applied_field_set){
         const double Hx=sim::H_vec[0]*sim::H_applied;
         const double Hy=sim::H_vec[1]*sim::H_applied;
         const double Hz=sim::H_vec[2]*sim::H_applied;
         for(int i=start_index;i<end_index; i++){
            fields_array_x[i]+=Hx;
            fields_array_y[i]+=Hy;
            fields_array_z[i]+=Hz;
         }
      }
   }

   // Collect all field terms of atom i that depends on s_i itself, once per single-spin update.
   // During the single-spin update, the position and neighbour spins are fixed so we gather pair
   // data once and reuse it for every rotation of the update. The pair data required is:
   // full Neel q:  bond e_ij, s_j and q(r_ij) for r_ij < r_cut_neel_q
   // biquadratic:  s_j and K(r_ij)/mu_i       for r_ij < r_cut_exchange
   // The anisotropy terms need no neighbours and are evaluated directly by anisotropy::nonlinear_spin_field()
   template<bool fully_periodic>
   void gather_nonlinear_field_terms(const int i,
                                     nonlinear_field_terms_t& terms){

      const unsigned int imat = atoms::type_array[i];
      terms.q_neighbours.clear();
      terms.biquadratic_neighbours.clear();
      terms.inverse_moment = 1.0/::mp::material[imat].mu_s_SI;
      terms.material = imat;

      const bool gather_q = sld::internal::full_neel_nonlinear_q;
      const bool gather_biquadratic = sld::internal::spin_hamiltonian == sld::internal::biquadratic_spin_hamiltonian;
      if(!gather_q && !gather_biquadratic) return; // nothing to gather

      // full Neel q settings
      const bool use_bethe_slater = sld::internal::neel_radial_function == sld::internal::bethe_slater_neel_radial_function;
      const bool use_smooth_cutoff = sld::internal::neel_cutoff_function == sld::internal::smooth_neel_cutoff_function;
      const double cutoff = sld::internal::r_cut_neel_q;
      const double cutoff_squared = cutoff*cutoff;

      // biquadratic exchange settings, as in compute_exchange()
      const bool use_bethe_slater_exchange = sld::internal::exchange_function == sld::internal::bethe_slater_exchange_function;
      const double exchange_cutoff_squared = sld::internal::r_cut_exchange*sld::internal::r_cut_exchange;
      const double exchange_inverse_cutoff = 1.0/sld::internal::r_cut_exchange;
      const double exchange_K0 = (gather_biquadratic && !use_bethe_slater_exchange) ? sld::internal::mp[imat].K0_ms.get() : 0.0;

      const int neighbour_start = atoms::neighbour_list_start_index[i];
      const int neighbour_end = atoms::neighbour_list_end_index[i]+1;
      for(int neighbour=neighbour_start; neighbour<neighbour_end; ++neighbour){
         const int j = atoms::neighbour_list_array[neighbour];
         if(j == i) continue;

         neel_vector_t displacement = {atoms::x_coord_array[i]-atoms::x_coord_array[j],
                                       atoms::y_coord_array[i]-atoms::y_coord_array[j],
                                       atoms::z_coord_array[i]-atoms::z_coord_array[j]};

         displacement.x = sld::PBC_wrap<fully_periodic>(displacement.x, cs::system_dimensions[0], cs::pbc[0]);
         displacement.y = sld::PBC_wrap<fully_periodic>(displacement.y, cs::system_dimensions[1], cs::pbc[1]);
         displacement.z = sld::PBC_wrap<fully_periodic>(displacement.z, cs::system_dimensions[2], cs::pbc[2]);

         const double distance_squared = neel_dot(displacement, displacement);
         const bool within_q = gather_q && distance_squared < cutoff_squared;
         const bool within_exchange = gather_biquadratic && distance_squared < exchange_cutoff_squared;
         if(!within_q && !within_exchange) continue;
         const double distance = std::sqrt(distance_squared);
         const neel_vector_t spin_j = {atoms::x_spin_array[j], atoms::y_spin_array[j], atoms::z_spin_array[j]};

         // biquadratic pair: K(r)(s_i.s_j)^2 gives the field 2(K/mu_i)(s_i.s_j)s_j
         if(within_exchange){
            double k;
            if(use_bethe_slater_exchange){
               k = bethe_slater(sld::internal::mp[imat].bethe_slater_alpha_k.get(),
                                sld::internal::mp[imat].bethe_slater_gamma_k.get(),
                                sld::internal::mp[imat].bethe_slater_delta_k.get(),
                                distance).value*terms.inverse_moment;
            }
            else{
               // K(r) = K0 (1 - r/r_c)^3, with K0_ms = K0/mu_i
               const double y = 1.0-distance*exchange_inverse_cutoff;
               k = exchange_K0*y*y*y;
            }
            terms.biquadratic_neighbours.push_back(biquadratic_neighbour_t{spin_j, k});
         }

         if(!within_q) continue;

         radial_result_t q_curve;
         if(use_bethe_slater){
            q_curve = bethe_slater(sld::internal::mp[imat].neel_alpha_q.get(),
                                   sld::internal::mp[imat].neel_gamma_q.get(),
                                   sld::internal::mp[imat].neel_delta_q.get(), distance);
         }
         else{
            q_curve = inverse_fourth(sld::internal::mp[imat].neel_C_q.get()*sld::internal::mp[imat].J0.get(), distance);
         }

         const unsigned int jmat = atoms::type_array[j];
         if(jmat != imat){
            radial_result_t q_curve_j;
            if(use_bethe_slater){
               q_curve_j = bethe_slater(sld::internal::mp[jmat].neel_alpha_q.get(),
                                        sld::internal::mp[jmat].neel_gamma_q.get(),
                                        sld::internal::mp[jmat].neel_delta_q.get(), distance);
            }
            else{
               q_curve_j = inverse_fourth(sld::internal::mp[jmat].neel_C_q.get()*sld::internal::mp[jmat].J0.get(), distance);
            }
            q_curve.value = 0.5*(q_curve.value+q_curve_j.value);
            q_curve.derivative = 0.5*(q_curve.derivative+q_curve_j.derivative);
         }

         if(use_smooth_cutoff){
            q_curve = apply_smooth_cutoff(q_curve, distance, sld::internal::r_switch_neel_q, cutoff);
         }
         else{
            q_curve = apply_hard_cutoff(q_curve, distance, cutoff);
         }
         if(q_curve.value == 0.0) continue;

         const double inverse_distance = 1.0/distance;
         const neel_vector_t bond = neel_scale(displacement, inverse_distance);
         terms.q_neighbours.push_back(neel_q_neighbour_t{bond, spin_j, q_curve.value});
      }

      return;
   }

   // Select the PBC mode once to avoid repeated boundary checks in the inner loops
   void gather_nonlinear_field_terms(const int i,
                                     nonlinear_field_terms_t& terms){
      const bool fully_periodic_boundaries = cs::pbc[0] && cs::pbc[1] && cs::pbc[2];
      if(fully_periodic_boundaries){
         gather_nonlinear_field_terms<true>(i, terms);
      }
      else{
         gather_nonlinear_field_terms<false>(i, terms);
      }
   }

   // Calculate only the q part of atom i's full-Neel spin field, including its
   // term linear in s_i, for a trial value s of that spin according to
   // B^q_i(s) = (1/mu_i) sum_j dPhi_q,ij/ds_i evaluated at s_i = s
   // sum in the same order as the gather above
   neel_vector_t compute_full_neel_q_field(const std::vector<neel_q_neighbour_t>& neighbours,
                                           const double inverse_moment,
                                           const neel_vector_t& trial_spin){

      neel_vector_t field = {0.0, 0.0, 0.0};
      for(size_t n=0; n<neighbours.size(); ++n){
         const neel_vector_t derivative = evaluate_neel_q_spin_derivative(neighbours[n].bond, trial_spin, neighbours[n].spin_j, neighbours[n].q);
         field = neel_add(field, neel_scale(derivative, inverse_moment));
      }
      return field;
   }

   // Calculate the part of atom i's field that depends on s_i itself, for a trial value s of that spin according to
   // B^nl_i(s) = B^q_i(s) + sum_j 2 (K_ij/mu_i)(s.s_j) s_j + B^an_i(s)
   // where B^an_i is the second order uniaxial plus fourth order cubic anisotropy field
   neel_vector_t compute_nonlinear_field(const nonlinear_field_terms_t& terms,
                                         const neel_vector_t& trial_spin){

      neel_vector_t field = compute_full_neel_q_field(terms.q_neighbours, terms.inverse_moment, trial_spin); // B^q_i(s)

      if(!terms.biquadratic_neighbours.empty()){
         neel_vector_t biquadratic = {0.0, 0.0, 0.0};
         for(size_t n=0; n<terms.biquadratic_neighbours.size(); ++n){
            const neel_vector_t& spin_j = terms.biquadratic_neighbours[n].spin_j;
            const double scale = 2.0*terms.biquadratic_neighbours[n].k*neel_dot(trial_spin, spin_j); // 2(K_ij/mu_i)(s.s_j)
            biquadratic = neel_add(biquadratic, neel_scale(spin_j, scale));
         }
         field = neel_add(field, biquadratic);
      }

      if(anisotropy::nonlinear_spin_field_enabled()){
         neel_vector_t anisotropy_field;
         anisotropy::nonlinear_spin_field(terms.material, trial_spin.x, trial_spin.y, trial_spin.z,
                                          anisotropy_field.x, anisotropy_field.y, anisotropy_field.z); // B^an_i(s)
         field = neel_add(field, anisotropy_field);
      }

      return field;
   }

   // Split B^nl_i(s) into fixed-axis terms g(a.s) a.
   // Each term is an energy F(a.s) of a single projection, so a rotation about a leaves a.s and the field constant.
   // This extends the single-spin decomposition of Weissenhofer et al., Phys. Rev. B 108, L060404 (2023) to the other non-linear terms supported.
   // Full Neel q: with x = e.s, y = e.s_j (fixed) and z = s_j.s, the pair polynomial
   //    Phi_q/q = (12/35)(xy - z/3) + (9/5)(x^2 - z/3)(y^2 - z/3) - (2/5)(xy^3 + yx^3)
   // is rewritten using x^2 z = [(x+z)^3 - (x-z)^3 - 2z^3]/6 as
   //    Phi_q/q = P_e(x) + P_j(z) - (x+z)^3/10 + (x-z)^3/10
   //    P_e(x) = (12y/35 - 2y^3/5) x + (9/5) y^2 x^2 - (2/5) y x^3
   //    P_j(z) = -(4/35 + 3y^2/5) z + z^2/5 + z^3/5
   // so with Q = q/mu_i the field Q grad(Phi_q/q) has the four axes e, s_j, e+s_j and e-s_j, each with g(u) = Q dP/du.
   // Biquadratic, uniaxial and cubic fields are already of this form, about s_j, the easy axis and the three Cartesian axes.
   void build_spin_rotation_terms(const nonlinear_field_terms_t& terms,
                                  std::vector<spin_rotation_term_t>& rotation_terms){

      rotation_terms.clear();

      // full Neel q, four terms per neighbour (axes e+s_j and e-s_j are not unit vectors)
      for(size_t n=0; n<terms.q_neighbours.size(); ++n){
         const neel_vector_t& e = terms.q_neighbours[n].bond;
         const neel_vector_t& t = terms.q_neighbours[n].spin_j;
         const double y = neel_dot(e, t);
         const double Q = terms.q_neighbours[n].q*terms.inverse_moment; // T
         rotation_terms.push_back(spin_rotation_term_t{e, Q*(12.0*y/35.0 - 0.4*y*y*y), 3.6*Q*y*y, -1.2*Q*y, 0.0});
         rotation_terms.push_back(spin_rotation_term_t{t, -Q*(4.0/35.0 + 0.6*y*y), 0.4*Q, 0.6*Q, 0.0});
         rotation_terms.push_back(spin_rotation_term_t{neel_add(e, t), 0.0, 0.0, -0.3*Q, 0.0});
         rotation_terms.push_back(spin_rotation_term_t{neel_add(e, neel_scale(t, -1.0)), 0.0, 0.0, 0.3*Q, 0.0});
      }

      // biquadratic exchange, field 2 (K_ij/mu_i)(s.s_j) s_j
      for(size_t n=0; n<terms.biquadratic_neighbours.size(); ++n){
         rotation_terms.push_back(spin_rotation_term_t{terms.biquadratic_neighbours[n].spin_j, 0.0, 2.0*terms.biquadratic_neighbours[n].k, 0.0, 0.0});
      }

      // second order uniaxial anisotropy, field 2 (K_u/mu_i)(s.e) e
      neel_vector_t easy_axis;
      double ku;
      if(anisotropy::uniaxial_second_order_coefficients(terms.material, easy_axis.x, easy_axis.y, easy_axis.z, ku)){
         rotation_terms.push_back(spin_rotation_term_t{easy_axis, 0.0, 2.0*ku, 0.0, 0.0});
      }

      // fourth order cubic anisotropy, field 2 (K_c/mu_i) s_a^3 along each Cartesian axis a
      double kc;
      if(anisotropy::cubic_fourth_order_coefficient(terms.material, kc)){
         rotation_terms.push_back(spin_rotation_term_t{neel_vector_t{1.0, 0.0, 0.0}, 0.0, 0.0, 0.0, 2.0*kc});
         rotation_terms.push_back(spin_rotation_term_t{neel_vector_t{0.0, 1.0, 0.0}, 0.0, 0.0, 0.0, 2.0*kc});
         rotation_terms.push_back(spin_rotation_term_t{neel_vector_t{0.0, 0.0, 1.0}, 0.0, 0.0, 0.0, 2.0*kc});
      }

      for(size_t n=0; n<rotation_terms.size(); ++n) rotation_terms[n].axis_squared = neel_dot(rotation_terms[n].axis, rotation_terms[n].axis);

      return;
   }

   // calculate distance dependent bilinear or biquadratic exchange fields and forces
   // spin_field_only stores only the field, so unused force/energy/curvature calculations are not done in this mode
   template<bool fully_periodic, bool spin_field_only>
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
               std::vector<double>& fields_array_x, //  vectors for fields
               std::vector<double>& fields_array_y,
               std::vector<double>& fields_array_z,
               nonlinear_field_terms_t& terms){ // filled only by the spin_field_only mode

       double rx, ry, rz;
       double dx, dy, dz;
       double sx, sy, sz;
       double sjx, sjy,sjz;
       double si_dot_sj;
       double fx = 0.0, fy = 0.0, fz = 0.0;
       double hx = 0.0, hy = 0.0, hz = 0.0;
       double rji_sqr, rji, inv_rji; //,  inv_rji2;
       double y, f_exch, energy = 0.0;
       //double exch_J0 = sld::internal::mp[0].J0_ms.get(); //7034.8836847351113; //
       //double exch_J0_prime = sld::internal::mp[0].J0_prime.get()/1.602176634e-19; //in J  0.72320000000000007 ;
       double J;
       int j;
       double r_sqr_cut=sld::internal::r_cut_exchange*sld::internal::r_cut_exchange;
       //double oneover3=1.0/3.0;
       double exch_inv_rcut=1.0/sld::internal::r_cut_exchange;
       double sumJ=0.0;
       double hessian_trace=0.0;
       const double joules_per_electron_volt = 1.602176634e-19;
       const bool use_bethe_slater =
          sld::internal::exchange_function ==
          sld::internal::bethe_slater_exchange_function;
       const bool use_biquadratic =
          sld::internal::spin_hamiltonian ==
          sld::internal::biquadratic_spin_hamiltonian;
       const double hamiltonian_offset =
          sld::internal::exchange_offset ? 1.0 : 0.0;
       // collect s_j and K(r)/mu_i of each biquadratic pair for the nonlinear hamiltonian single spin update
       const bool collect_biquadratic = spin_field_only && use_biquadratic && sld::internal::nonlinear_spin_hamiltonian;

       for(int i=start_index;i<end_index; ++i){

          const unsigned int imat = atoms::type_array[i];
          double exch_J0 = sld::internal::mp[imat].J0_ms.get(); //7034.8836847351113; //
          double exch_J0_prime = sld::internal::mp[imat].J0_prime.get()/joules_per_electron_volt;
          double exch_K0 = 0.0;
          double exch_K0_prime = 0.0;
          if(use_biquadratic && !use_bethe_slater){
             exch_K0 = sld::internal::mp[imat].K0_ms.get();
             exch_K0_prime =
                sld::internal::mp[imat].K0_prime.get()/joules_per_electron_volt;
          }
          const double inverse_moment = 1.0/::mp::material[imat].mu_s_SI;
          //int count_int=0;

          fx = 0.0;
          fy = 0.0;
          fz = 0.0;
          hx=0.0;
          hy=0.0;
          hz=0.0;
          sumJ=0.0;
          hessian_trace=0.0;
          energy=0.0;

          rx = x_coord_array[i];
          ry = y_coord_array[i];
          rz = z_coord_array[i];

          sx = x_spin_array[i];
          sy = y_spin_array[i];
          sz = z_spin_array[i];

          int nbr_start = neighbour_list_start_index[i];
          int nbr_end = neighbour_list_end_index[i]+1;

          for( int n = nbr_start; n < nbr_end; ++n){

            j = neighbour_list_array[n];

            if ( j != i){
             dx = -x_coord_array[j] + rx;
             dy = -y_coord_array[j] + ry;
             dz = -z_coord_array[j] + rz;

             dx = sld::PBC_wrap<fully_periodic>(dx, cs::system_dimensions[0], cs::pbc[0]);
             dy = sld::PBC_wrap<fully_periodic>(dy, cs::system_dimensions[1], cs::pbc[1]);
             dz = sld::PBC_wrap<fully_periodic>(dz, cs::system_dimensions[2], cs::pbc[2]);


             rji_sqr = dx*dx + dy*dy + dz*dz;

             if( rji_sqr < r_sqr_cut)
             {   //count_int++;

                 rji = sqrt(rji_sqr);
                 inv_rji = 1.0/ rji;

                 sjx = x_spin_array[j];
                 sjy = y_spin_array[j];
                 sjz = z_spin_array[j];

                 //this part calculates the exchange forces
                 //for computational efficiency, forces and fields are calculated at the same time
                 si_dot_sj = sx * sjx + sy * sjy + sz * sjz;

                 if(use_bethe_slater){
                    const radial_result_t j_curve =
                       bethe_slater(sld::internal::mp[imat].bethe_slater_alpha_j.get(),
                                    sld::internal::mp[imat].bethe_slater_gamma_j.get(),
                                    sld::internal::mp[imat].bethe_slater_delta_j.get(),
                                    rji);
                    const double j_field = j_curve.value*inverse_moment;
                    double pair_field = j_field;
                    double radial_derivative = j_curve.derivative*(si_dot_sj-hamiltonian_offset);
                    energy += j_field*(si_dot_sj-hamiltonian_offset);

                    // The biquadratic mode adds K(r)(si.sj)^2 to the Bethe-Slater J(r) Hamiltonian.
                    if(use_biquadratic){
                       const radial_result_t k_curve =
                          bethe_slater(sld::internal::mp[imat].bethe_slater_alpha_k.get(),
                                       sld::internal::mp[imat].bethe_slater_gamma_k.get(),
                                       sld::internal::mp[imat].bethe_slater_delta_k.get(),
                                       rji);
                       const double k_field = k_curve.value*inverse_moment;
                       if(collect_biquadratic) terms.biquadratic_neighbours.push_back(biquadratic_neighbour_t{neel_vector_t{sjx, sjy, sjz}, k_field});
                       pair_field += 2.0*k_field*si_dot_sj;
                       // Tr[(I-s_i s_i^T)dB_i^K/ds_i]=2(K/mu_i)[1-(s_i.s_j)^2]
                       hessian_trace +=
                          2.0*k_field*(1.0-si_dot_sj*si_dot_sj);
                       radial_derivative +=
                          k_curve.derivative*
                          (si_dot_sj*si_dot_sj-hamiltonian_offset);
                       energy +=
                          k_field*(si_dot_sj*si_dot_sj-hamiltonian_offset);
                    }

                    hx += pair_field*sjx;
                    hy += pair_field*sjy;
                    hz += pair_field*sjz;
                    sumJ += pair_field;

                    // The optional offsets remove collinear ground state forces without changing spin precession
                    f_exch = radial_derivative/joules_per_electron_volt;
                    fx += f_exch*dx*inv_rji;
                    fy += f_exch*dy*inv_rji;
                    fz += f_exch*dz*inv_rji;
                 }
                 else{
                    y = (1.0 - rji * exch_inv_rcut);
                    J = exch_J0 * y * y * y;
                    double pair_field = J;
                    double radial_derivative =
                       -exch_J0_prime*y*y*
                       (si_dot_sj-hamiltonian_offset);
                    energy += J*(si_dot_sj-hamiltonian_offset);

                    if(use_biquadratic){
                       const double K = exch_K0*y*y*y;
                       if(collect_biquadratic) terms.biquadratic_neighbours.push_back(biquadratic_neighbour_t{neel_vector_t{sjx, sjy, sjz}, K});
                       pair_field += 2.0*K*si_dot_sj;
                       hessian_trace += 2.0*K*(1.0-si_dot_sj*si_dot_sj);
                       radial_derivative +=
                          -exch_K0_prime*y*y*
                          (si_dot_sj*si_dot_sj-hamiltonian_offset);
                       energy +=
                          K*(si_dot_sj*si_dot_sj-hamiltonian_offset);
                    }

                    hx += pair_field*sjx;
                    hy += pair_field*sjy;
                    hz += pair_field*sjz;
                    sumJ += pair_field;

                    f_exch = radial_derivative;
                    fx += f_exch*dx*inv_rji;
                    fy += f_exch*dy*inv_rji;
                    fz += f_exch*dz*inv_rji;
                 }
                 //std::cout<<std::setprecision(15)<<std::endl;
                 /*
                 if(abs(x_coord_array[i]-20.09)<1e-3 &&abs(y_coord_array[i]-20.09)<1e-3  && abs(z_coord_array[i]-20.09)<1e-3 ){
                  std::cout<<"i= "<<i<<" j= "<<j<<" pos  "<<"\t" << rx<<"\t"<<ry<<"\t"<<rz<<"\t"<<x_coord_array[j]<<"\t"<<y_coord_array[j]<<"\t"<<z_coord_array[j]<<std::endl;
                 std::cout<<"i= "<<i<<" j= "<<j<<" spin  "<<"\t" << sx<<"\t"<<sy<<"\t"<<sz<<"\t"<<sjx<<"\t"<<sjy<<"\t"<<sjz<<std::endl;
                 std::cout<<i<<" Exchange:   Forces: "<<f_exch * dx *  (si_dot_sj)* inv_rji<<"\t"<<f_exch * dy *  (si_dot_sj)* inv_rji<<"\t"<<f_exch * dz *  (si_dot_sj)* inv_rji<<std::endl;
                 std::cout<<i<<" Exchange:   Fields: "<<J * sjx<<"\t"<<J * sjy<<"\t"<<J * sjz<<std::endl;
                 std::cout<<i<<" Exchange:   dxyz: "<<dx<<"\t"<<dy<<"\t"<<dz<<std::endl;

                 //std::cout<<"i= "<<i<<" j= "<<j<<"\t"<<" y= "<<y<<"\t"<<"exch_J0_prime: "<<exch_J0_prime<<"exch_J0 "<<exch_J0<<std::endl;
                 //std::cout<<i<<"energ "<<"\t"<<energy<<std::endl;
                 }*/

             }//end if cutoff

          }//end if i!=j

       }//end for loop neighbourlist

       if(!spin_field_only){
          forces_array_x[i] += fx;
          forces_array_y[i] += fy;
          forces_array_z[i] += fz;
       }

       fields_array_x[i] += hx;
       fields_array_y[i] += hy;
       fields_array_z[i] += hz;

       /*
       std::cout<<"i= "<<i<<" pos  "<<"\t" << rx<<"\t"<<ry<<"\t"<<rz<<std::endl;
       std::cout<<i<<"energ "<<"\t"<<energy<<std::endl;
       std::cout<<i<<"\t"<<hx<<"\t"<<hy<<"\t"<<hz<<"\t"<<fx<<"\t"<<fy<<"\t"<<fz<<std::endl;*/
       if(!spin_field_only){
          sld::internal::sumJ[i]=sumJ;
          sld::internal::exch_eng[i]=-0.5*energy;
          sld::internal::spin_hessian_trace[i] += hessian_trace;
       }
       //std::cout<<"exchange int"<<count_int<<std::endl;

       /*if(abs(x_coord_array[i]-20.09)<1e-3 &&abs(y_coord_array[i]-20.09)<1e-3  && abs(z_coord_array[i]-20.09)<1e-3 ){
       std::cout<<std::setprecision(15)<<std::endl;
       std::cout<<"forces exch "<<forces_array_x[i]<<"\t"<<forces_array_y[i]<<"\t"<<forces_array_z[i]<<"\t"<<std::endl;
                    std::cout<<"xyz "<<x_coord_array[i]<<"\t"<<y_coord_array[i]<<"\t"<<z_coord_array[i]<<std::endl;
                    std::cout<<"sxyz "<<x_spin_array[i]<<"\t"<<y_spin_array[i]<<"\t"<<z_spin_array[i]<<std::endl;


         }*/


   }
   return;


}

   void compute_exchange(const int start_index,
               const int end_index,
               const std::vector<int>& neighbour_list_start_index,
               const std::vector<int>& neighbour_list_end_index,
               const std::vector<int>& type_array,
               const std::vector<int>& neighbour_list_array,
               const std::vector<double>& x_coord_array,
               const std::vector<double>& y_coord_array,
               const std::vector<double>& z_coord_array,
               const std::vector<double>& x_spin_array,
               const std::vector<double>& y_spin_array,
               const std::vector<double>& z_spin_array,
               std::vector<double>& forces_array_x,
               std::vector<double>& forces_array_y,
               std::vector<double>& forces_array_z,
               std::vector<double>& fields_array_x,
               std::vector<double>& fields_array_y,
               std::vector<double>& fields_array_z){

      // Select the PBC mode once to avoid repeated boundary checks in the inner loops
      const bool fully_periodic_boundaries = cs::pbc[0] && cs::pbc[1] && cs::pbc[2];
      nonlinear_field_terms_t unused_terms; // nonlinear terms are collected by compute_spin_field
      if(fully_periodic_boundaries){
         compute_exchange<true, false>(
            start_index, end_index, neighbour_list_start_index,
            neighbour_list_end_index, type_array, neighbour_list_array,
            x_coord_array, y_coord_array, z_coord_array,
            x_spin_array, y_spin_array, z_spin_array,
            forces_array_x, forces_array_y, forces_array_z,
            fields_array_x, fields_array_y, fields_array_z, unused_terms);
      }
      else{
         compute_exchange<false, false>(
            start_index, end_index, neighbour_list_start_index,
            neighbour_list_end_index, type_array, neighbour_list_array,
            x_coord_array, y_coord_array, z_coord_array,
            x_spin_array, y_spin_array, z_spin_array,
            forces_array_x, forces_array_y, forces_array_z,
            fields_array_x, fields_array_y, fields_array_z, unused_terms);
      }
   }

   // calculate the pseudodipolar spin-lattice coupling fields and forces
   // spin_field_only stores only the field and coupling field to reduce uneccessary calculations 
   template<bool fully_periodic, bool spin_field_only>
   void compute_sld_coupling(const int start_index,
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
               std::vector<double>& fields_array_x, //  vectors for fields
               std::vector<double>& fields_array_y,
               std::vector<double>& fields_array_z){


                  double rx, ry, rz;
                  double dx, dy, dz;
                  double sx, sy, sz;
                  double sjx, sjy,sjz;
                  double si_dot_sj;
                  double fc_x = 0.0, fc_y = 0.0, fc_z = 0.0;
                  double hc_x = 0.0, hc_y = 0.0, hc_z = 0.0;
                  double rji_sqr, rji, inv_rji,  inv_rji2, inv_rji4, inv_rji6;
                  double sj_dot_rji, si_dot_rji;
                  double energy_c;
                  int j; //, count_int;

                  double r_sqr_cut=sld::internal::r_cut_fields*sld::internal::r_cut_fields;
                  double oneover3=1.0/3.0;
                  double sumC;

                  for(int i=start_index;i<end_index; ++i){

                     const unsigned int imat = atoms::type_array[i];
                     double fact=sld::internal::mp[imat].C0.get()/1.602176634e-19;//in J, 0.4520;//
                     double fact_ms= sld::internal::mp[imat].C0_ms.get();//3517.4418423675556;

                     //count_int=0;
                     fc_x = 0.0;
                     fc_y = 0.0;
                     fc_z = 0.0;
                     hc_x = 0.0;
                     hc_y = 0.0;
                     hc_z = 0.0;
                     energy_c=0.0;
                     sumC=0.0;

                     rx = x_coord_array[i];
                     ry = y_coord_array[i];
                     rz = z_coord_array[i];

                     sx = x_spin_array[i];
                     sy = y_spin_array[i];
                     sz = z_spin_array[i];

                     int nbr_start = neighbour_list_start_index[i];
                     int nbr_end = neighbour_list_end_index[i]+1;


                     for( int n = nbr_start; n < nbr_end; ++n){
                       j = neighbour_list_array[n];


                       if ( j != i){

                       dx = -x_coord_array[j] + rx;
                       dy = -y_coord_array[j] + ry;
                       dz = -z_coord_array[j] + rz;

                       dx = sld::PBC_wrap<fully_periodic>(dx, cs::system_dimensions[0], cs::pbc[0]);
                       dy = sld::PBC_wrap<fully_periodic>(dy, cs::system_dimensions[1], cs::pbc[1]);
                       dz = sld::PBC_wrap<fully_periodic>(dz, cs::system_dimensions[2], cs::pbc[2]);


                       rji_sqr = dx*dx + dy*dy + dz*dz;

                       if( rji_sqr < r_sqr_cut)
                       {

                            //count_int++;
                            rji = sqrt(rji_sqr);
                            inv_rji = 1.0/ rji;


                            sjx = x_spin_array[j];
                            sjy = y_spin_array[j];
                            sjz = z_spin_array[j];


                            //std::cout<<"spin "<<sx<<"\t"<<sy<<"\t"<<sz<<"\t"<<sjx<<"\t"<<sjy<<"\t"<<sjz<<std::endl;
                            si_dot_sj = sx * sjx + sy * sjy + sz * sjz;

                            sj_dot_rji = (dx * sjx + dy * sjy + dz * sjz);
                            si_dot_rji = (dx * sx  + dy * sy  + dz * sz);

                            inv_rji2=inv_rji*inv_rji;
                            inv_rji4=inv_rji2*inv_rji2;
                            inv_rji6=inv_rji4*inv_rji2;


                            //adding fields from the pseudo-dipolar coupling
                            hc_x +=fact_ms*inv_rji4*(inv_rji2*dx*sj_dot_rji-oneover3*sjx);
                            hc_y +=fact_ms*inv_rji4*(inv_rji2*dy*sj_dot_rji-oneover3*sjy);
                            hc_z +=fact_ms*inv_rji4*(inv_rji2*dz*sj_dot_rji-oneover3*sjz);

                            // phi_pseudodp=f(r)[(e.s_i)(e.s_j)-(s_i.s_j)/3], factor 1/2 to remove double counting from the neighbour list
                            energy_c += fact_ms*inv_rji4*(inv_rji2*si_dot_rji*sj_dot_rji-oneover3*si_dot_sj);

                            //adding forces from the pseudo-dipolar coupling
                            fc_x += fact*inv_rji6*( sj_dot_rji * sx + si_dot_rji * sjx -6.0* dx* sj_dot_rji* si_dot_rji * inv_rji2+ oneover3*4.0*si_dot_sj*dx);
                            fc_y += fact*inv_rji6*( sj_dot_rji * sy + si_dot_rji * sjy -6.0* dy* sj_dot_rji* si_dot_rji * inv_rji2+ oneover3*4.0*si_dot_sj*dy);
                            fc_z += fact*inv_rji6*( sj_dot_rji * sz + si_dot_rji * sjz -6.0* dz* sj_dot_rji* si_dot_rji * inv_rji2+ oneover3*4.0*si_dot_sj*dz);



                            sumC +=fact_ms*inv_rji4;

                           /*
                            double fc_x1 = fact*inv_rji6*( sj_dot_rji * sx + si_dot_rji * sjx -6.0* dx* sj_dot_rji* si_dot_rji * inv_rji2+ oneover3*4.0*si_dot_sj*dx);
                            double fc_y1 = fact*inv_rji6*( sj_dot_rji * sy + si_dot_rji * sjy -6.0* dy* sj_dot_rji* si_dot_rji * inv_rji2+ oneover3*4.0*si_dot_sj*dy);
                            double fc_z1 = fact*inv_rji6*( sj_dot_rji * sz + si_dot_rji * sjz -6.0* dz* sj_dot_rji* si_dot_rji * inv_rji2+ oneover3*4.0*si_dot_sj*dz);

                            //if(i==1110) ofile_cp<<std::setprecision(17)<< rx<<"\t"<<ry<<"\t"<<rz<<"\t"<<i<<"\t" <<j<<"\t"<<"\t"<<x_coord_array[j]<<"\t"<<y_coord_array[j]<<"\t"<<z_coord_array[j]<<"\t"<< hc_x1 << "\t" << hc_y1<< "\t" <<hc_z1<<"\t"<< fc_x1 << "\t" << fc_y1<< "\t" <<fc_z1<<std::endl;

                           std::cout<<"i="<<i<<" j="<<j<<std::endl;
                           std::cout<<i<<" rji "<<"\t"<<rji<<"\t"<<inv_rji<<std::endl;
                           std::cout<<i<<" pos  "<<"\t" << rx<<"\t"<<ry<<"\t"<<rz<<"\t"<<x_coord_array[j]<<"\t"<<y_coord_array[j]<<"\t"<<z_coord_array[j]<<std::endl;
                           std::cout<<i<<" forces coup " << fc_x1 << "\t" << fc_y1<< "\t" <<fc_z1<<std::endl;
                           std::cout<<i<<" fields coup " << hc_x1 << "\t" << hc_y1<< "\t" <<hc_z1<<std::endl;
                           std::cout<<i<<"energ "<<"\t"<<energy_c<<"\t"<<count_int<<std::endl;*/
                  }
               }


            }



            if(!spin_field_only){
               forces_array_x[i] += fc_x;
               forces_array_y[i] += fc_y;
               forces_array_z[i] += fc_z;
            }

            fields_array_x[i] += hc_x;
            fields_array_y[i] += hc_y;
            fields_array_z[i] += hc_z;

            sld::internal::coupling_field_x[i] = hc_x;
            sld::internal::coupling_field_y[i] = hc_y;
            sld::internal::coupling_field_z[i] = hc_z;

            if(!spin_field_only){
               sld::internal::sumC[i]=sumC;
               sld::internal::coupl_eng[i]=-0.5*energy_c;
            }

            }


            return;
         } // end function compute_sld_coupling

   void compute_sld_coupling(const int start_index,
               const int end_index,
               const std::vector<int>& neighbour_list_start_index,
               const std::vector<int>& neighbour_list_end_index,
               const std::vector<int>& type_array,
               const std::vector<int>& neighbour_list_array,
               const std::vector<double>& x_coord_array,
               const std::vector<double>& y_coord_array,
               const std::vector<double>& z_coord_array,
               const std::vector<double>& x_spin_array,
               const std::vector<double>& y_spin_array,
               const std::vector<double>& z_spin_array,
               std::vector<double>& forces_array_x,
               std::vector<double>& forces_array_y,
               std::vector<double>& forces_array_z,
               std::vector<double>& fields_array_x,
               std::vector<double>& fields_array_y,
               std::vector<double>& fields_array_z){

      // Select the PBC mode once to avoid repeated boundary checks in the inner loops
      const bool fully_periodic_boundaries = cs::pbc[0] && cs::pbc[1] && cs::pbc[2];
      if(fully_periodic_boundaries){
         compute_sld_coupling<true, false>(
            start_index, end_index, neighbour_list_start_index,
            neighbour_list_end_index, type_array, neighbour_list_array,
            x_coord_array, y_coord_array, z_coord_array,
            x_spin_array, y_spin_array, z_spin_array,
            forces_array_x, forces_array_y, forces_array_z,
            fields_array_x, fields_array_y, fields_array_z);
      }
      else{
         compute_sld_coupling<false, false>(
            start_index, end_index, neighbour_list_start_index,
            neighbour_list_end_index, type_array, neighbour_list_array,
            x_coord_array, y_coord_array, z_coord_array,
            x_spin_array, y_spin_array, z_spin_array,
            forces_array_x, forces_array_y, forces_array_z,
            fields_array_x, fields_array_y, fields_array_z);
      }
   }


   // calculate full Neel coupling fields, forces, energies and spin curvature
   // spin_field_only evaluates only the spin derivative of each pair and stores only the field and coupling field
   template<bool fully_periodic, bool spin_field_only>
   void compute_sld_coupling_neel(const int start_index,
               const int end_index, // last +1 atom to be calculated
               const std::vector<int>& neighbour_list_start_index,
               const std::vector<int>& neighbour_list_end_index,
               const std::vector<int>& type_array, // type for atom
               const std::vector<int>& neighbour_list_array, // list of interactions between atom
               const std::vector<double>& x_coord_array, // coord vectors for atoms
               const std::vector<double>& y_coord_array,
               const std::vector<double>& z_coord_array,
               const std::vector<double>& x_spin_array, // coord vectors for spins
               const std::vector<double>& y_spin_array,
               const std::vector<double>& z_spin_array,
               std::vector<double>& forces_array_x, //  vectors for forces
               std::vector<double>& forces_array_y,
               std::vector<double>& forces_array_z,
               std::vector<double>& fields_array_x, //  vectors for fields
               std::vector<double>& fields_array_y,
               std::vector<double>& fields_array_z,
               nonlinear_field_terms_t& terms){ // filled only by the spin_field_only mode

      const double joules_per_electron_volt = 1.602176634e-19;
      const double maximum_cutoff =
         std::max(sld::internal::r_cut_neel_l, sld::internal::r_cut_neel_q);
      const double maximum_cutoff_squared = maximum_cutoff*maximum_cutoff;
      const bool use_bethe_slater =
         sld::internal::neel_radial_function ==
         sld::internal::bethe_slater_neel_radial_function;
      const bool use_smooth_cutoff =
         sld::internal::neel_cutoff_function ==
         sld::internal::smooth_neel_cutoff_function;
      // collect bond, s_j and q(r) of each q pair for the nonlinear single spin update. Follows gather_nonlinear_field_terms where a pair is
      // added if r^2 < r_cut_q^2 and a nonzero q after the cutoff, so neighbours within only the l cutoff add no q term
      const bool collect_q = spin_field_only && sld::internal::full_neel_nonlinear_q;
      const double q_cutoff_squared = sld::internal::r_cut_neel_q*sld::internal::r_cut_neel_q;

      // H_Neel=-1/2 sum_ij Phi_ij, Phi_ij=l1*A+q1*B+q2*C
      for(int i=start_index; i<end_index; ++i){

         const unsigned int imat = type_array[i];
         const double inverse_moment = 1.0/::mp::material[imat].mu_s_SI;
         const neel_vector_t spin_i = {
            x_spin_array[i], y_spin_array[i], z_spin_array[i]
         };
         double force_x = 0.0;
         double force_y = 0.0;
         double force_z = 0.0;
         double field_x = 0.0;
         double field_y = 0.0;
         double field_z = 0.0;
         double hessian_trace = 0.0;
         double ordered_energy = 0.0;

         const int neighbour_start = neighbour_list_start_index[i];
         const int neighbour_end = neighbour_list_end_index[i]+1;
         for(int neighbour=neighbour_start; neighbour<neighbour_end; ++neighbour){

            const int j = neighbour_list_array[neighbour];
            if(j == i) continue;

            neel_vector_t displacement = {
               x_coord_array[i]-x_coord_array[j],
               y_coord_array[i]-y_coord_array[j],
               z_coord_array[i]-z_coord_array[j]
            };
            displacement.x = sld::PBC_wrap<fully_periodic>(displacement.x, cs::system_dimensions[0], cs::pbc[0]);
            displacement.y = sld::PBC_wrap<fully_periodic>(displacement.y, cs::system_dimensions[1], cs::pbc[1]);
            displacement.z = sld::PBC_wrap<fully_periodic>(displacement.z, cs::system_dimensions[2], cs::pbc[2]);

            const double distance_squared = neel_dot(displacement, displacement);
            if(distance_squared >= maximum_cutoff_squared) continue;
            const double distance = std::sqrt(distance_squared);

            radial_result_t l_curve;
            radial_result_t q_curve;
            // calculate Bethe-Slater or inverse fourth Neel radial functions
            if(use_bethe_slater){
               l_curve = bethe_slater(
                  sld::internal::mp[imat].neel_alpha_l.get(),
                  sld::internal::mp[imat].neel_gamma_l.get(),
                  sld::internal::mp[imat].neel_delta_l.get(),
                  distance);
               q_curve = bethe_slater(
                  sld::internal::mp[imat].neel_alpha_q.get(),
                  sld::internal::mp[imat].neel_gamma_q.get(),
                  sld::internal::mp[imat].neel_delta_q.get(),
                  distance);
            }
            else{
               // l=C_l*J0*(1 Angstrom/r)^4
               // q=C_q*J0*(1 Angstrom/r)^4
               l_curve = inverse_fourth(
                  sld::internal::mp[imat].neel_C_l.get()*
                  sld::internal::mp[imat].J0.get(), distance);
               q_curve = inverse_fourth(
                  sld::internal::mp[imat].neel_C_q.get()*
                  sld::internal::mp[imat].J0.get(), distance);
            }

            const unsigned int jmat = type_array[j];
            if(jmat != imat){
               radial_result_t l_curve_j;
               radial_result_t q_curve_j;
               if(use_bethe_slater){
                  l_curve_j = bethe_slater(
                     sld::internal::mp[jmat].neel_alpha_l.get(),
                     sld::internal::mp[jmat].neel_gamma_l.get(),
                     sld::internal::mp[jmat].neel_delta_l.get(),
                     distance);
                  q_curve_j = bethe_slater(
                     sld::internal::mp[jmat].neel_alpha_q.get(),
                     sld::internal::mp[jmat].neel_gamma_q.get(),
                     sld::internal::mp[jmat].neel_delta_q.get(),
                     distance);
               }
               else{
                  l_curve_j = inverse_fourth(
                     sld::internal::mp[jmat].neel_C_l.get()*
                     sld::internal::mp[jmat].J0.get(), distance);
                  q_curve_j = inverse_fourth(
                     sld::internal::mp[jmat].neel_C_q.get()*
                     sld::internal::mp[jmat].J0.get(), distance);
               }

               // symmetrise curves for different materials so energy and force remain conservative
               l_curve.value = 0.5*(l_curve.value+l_curve_j.value);
               l_curve.derivative = 0.5*(l_curve.derivative+l_curve_j.derivative);
               q_curve.value = 0.5*(q_curve.value+q_curve_j.value);
               q_curve.derivative = 0.5*(q_curve.derivative+q_curve_j.derivative);
            }

            // apply the selected hard or smooth cutoff independently to l and q
            if(use_smooth_cutoff){
               l_curve = apply_smooth_cutoff(
                  l_curve, distance, sld::internal::r_switch_neel_l,
                  sld::internal::r_cut_neel_l);
               q_curve = apply_smooth_cutoff(
                  q_curve, distance, sld::internal::r_switch_neel_q,
                  sld::internal::r_cut_neel_q);
            }
            else{
               l_curve = apply_hard_cutoff(
                  l_curve, distance, sld::internal::r_cut_neel_l);
               q_curve = apply_hard_cutoff(
                  q_curve, distance, sld::internal::r_cut_neel_q);
            }

            if(l_curve.value == 0.0 && l_curve.derivative == 0.0 &&
               q_curve.value == 0.0 && q_curve.derivative == 0.0) continue;

            const neel_vector_t spin_j = {
               x_spin_array[j], y_spin_array[j], z_spin_array[j]
            };
            const neel_pair_result_t pair =
               evaluate_neel_pair<spin_field_only>(displacement, spin_i, spin_j, l_curve, q_curve);
            if(collect_q && distance_squared < q_cutoff_squared && q_curve.value != 0.0){
               terms.q_neighbours.push_back(neel_q_neighbour_t{neel_scale(displacement, 1.0/distance), spin_j, q_curve.value});
            }

            // H_i=-(1/mu_i)dH/ds_i=(1/mu_i)sum_j dPhi_ij/ds_i
            field_x += pair.spin_derivative_i.x*inverse_moment;
            field_y += pair.spin_derivative_i.y*inverse_moment;
            field_z += pair.spin_derivative_i.z*inverse_moment;
            hessian_trace += pair.spin_hessian_trace_i*inverse_moment;

            // F_i=-dH/dr_i=sum_j grad_i(Phi_ij) and convert from J/A to eV/A
            force_x += pair.force_i.x/joules_per_electron_volt;
            force_y += pair.force_i.y/joules_per_electron_volt;
            force_z += pair.force_i.z/joules_per_electron_volt;
            ordered_energy += pair.energy;
         }

         // accumulate forces and fields for each atom
         if(!spin_field_only){
            forces_array_x[i] += force_x;
            forces_array_y[i] += force_y;
            forces_array_z[i] += force_z;
         }
         fields_array_x[i] += field_x;
         fields_array_y[i] += field_y;
         fields_array_z[i] += field_z;

         sld::internal::coupling_field_x[i] = field_x;
         sld::internal::coupling_field_y[i] = field_y;
         sld::internal::coupling_field_z[i] = field_z;
         if(!spin_field_only){
            sld::internal::spin_hessian_trace[i] += hessian_trace;
            sld::internal::sumC[i] = field_x*spin_i.x+field_y*spin_i.y+field_z*spin_i.z;

            // Each physical pair occurs twice, so store -1/2 sum_j Phi_ij in field units
            sld::internal::coupl_eng[i] = -0.5*ordered_energy*inverse_moment;
         }
      }

      return;
   } // end function compute_sld_coupling_neel

   void compute_sld_coupling_neel(const int start_index,
               const int end_index,
               const std::vector<int>& neighbour_list_start_index,
               const std::vector<int>& neighbour_list_end_index,
               const std::vector<int>& type_array,
               const std::vector<int>& neighbour_list_array,
               const std::vector<double>& x_coord_array,
               const std::vector<double>& y_coord_array,
               const std::vector<double>& z_coord_array,
               const std::vector<double>& x_spin_array,
               const std::vector<double>& y_spin_array,
               const std::vector<double>& z_spin_array,
               std::vector<double>& forces_array_x,
               std::vector<double>& forces_array_y,
               std::vector<double>& forces_array_z,
               std::vector<double>& fields_array_x,
               std::vector<double>& fields_array_y,
               std::vector<double>& fields_array_z){

      // Select the PBC mode once to avoid repeated boundary checks in the inner loops
      const bool fully_periodic_boundaries = cs::pbc[0] && cs::pbc[1] && cs::pbc[2];
      nonlinear_field_terms_t unused_terms;
      if(fully_periodic_boundaries){
         compute_sld_coupling_neel<true, false>(
            start_index, end_index, neighbour_list_start_index,
            neighbour_list_end_index, type_array, neighbour_list_array,
            x_coord_array, y_coord_array, z_coord_array,
            x_spin_array, y_spin_array, z_spin_array,
            forces_array_x, forces_array_y, forces_array_z,
            fields_array_x, fields_array_y, fields_array_z, unused_terms);
      }
      else{
         compute_sld_coupling_neel<false, false>(
            start_index, end_index, neighbour_list_start_index,
            neighbour_list_end_index, type_array, neighbour_list_array,
            x_coord_array, y_coord_array, z_coord_array,
            x_spin_array, y_spin_array, z_spin_array,
            forces_array_x, forces_array_y, forces_array_z,
            fields_array_x, fields_array_y, fields_array_z, unused_terms);
      }
   }

   // Field only version of compute_fields for one atom, mechanical forces, energies, sumJ, sumC and curvature are not calculated
   // because every function that uses these values first recalculates them with compute_fields
   // for nonlinear hamiltonians, exchange and full Neel loops loop the same neighbours so the neighbours are not traversed a second time
   template<bool fully_periodic>
   void compute_spin_field(const int atom, nonlinear_field_terms_t& terms){

      const unsigned int imat = atoms::type_array[atom];
      terms.q_neighbours.clear();
      terms.biquadratic_neighbours.clear();
      terms.inverse_moment = 1.0/::mp::material[imat].mu_s_SI;
      terms.material = imat;

      compute_exchange<fully_periodic, true>(atom, atom+1, atoms::neighbour_list_start_index, atoms::neighbour_list_end_index,
                                             atoms::type_array, atoms::neighbour_list_array,
                                             atoms::x_coord_array, atoms::y_coord_array, atoms::z_coord_array,
                                             atoms::x_spin_array, atoms::y_spin_array, atoms::z_spin_array,
                                             forces_array_x, forces_array_y, forces_array_z,
                                             fields_array_x, fields_array_y, fields_array_z, terms);

      if(pseudodipolar){ compute_sld_coupling<fully_periodic, true>(atom, atom+1, atoms::neighbour_list_start_index, atoms::neighbour_list_end_index,
                                                                   atoms::type_array, atoms::neighbour_list_array,
                                                                   atoms::x_coord_array, atoms::y_coord_array, atoms::z_coord_array,
                                                                   atoms::x_spin_array, atoms::y_spin_array, atoms::z_spin_array,
                                                                   forces_array_x, forces_array_y, forces_array_z,
                                                                   fields_array_x, fields_array_y, fields_array_z);
      }

      if(full_neel){ compute_sld_coupling_neel<fully_periodic, true>(atom, atom+1, atoms::neighbour_list_start_index, atoms::neighbour_list_end_index,
                                                                    atoms::type_array, atoms::neighbour_list_array,
                                                                    atoms::x_coord_array, atoms::y_coord_array, atoms::z_coord_array,
                                                                    atoms::x_spin_array, atoms::y_spin_array, atoms::z_spin_array,
                                                                    forces_array_x, forces_array_y, forces_array_z,
                                                                    fields_array_x, fields_array_y, fields_array_z, terms);
      }

      add_applied_field(atom, atom+1, fields_array_x, fields_array_y, fields_array_z);

      anisotropy::fields(atoms::x_spin_array, atoms::y_spin_array, atoms::z_spin_array, atoms::type_array,
                         fields_array_x, fields_array_y, fields_array_z, atom, atom+1, sim::temperature);
      return;
   }

   // Select the PBC mode once to avoid repeated boundary checks in the inner loops
   void compute_spin_field(const int atom, nonlinear_field_terms_t& terms){
      const bool fully_periodic_boundaries = cs::pbc[0] && cs::pbc[1] && cs::pbc[2];
      if(fully_periodic_boundaries){
         compute_spin_field<true>(atom, terms);
      }
      else{
         compute_spin_field<false>(atom, terms);
      }
   }

} // end of internal namespace
} // end of sld namespace
