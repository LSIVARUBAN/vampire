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
#include <iostream>
#include <cmath>
#include <vector>

// Vampire headers
#include "atoms.hpp"
#include "create.hpp"
#include "errors.hpp"
#include "material.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "sld.hpp"

//sld module headers M Strungaru
#include "internal.hpp"

namespace sld{
   void compute_forces_fields(const int start_index, // first atom for exchange interactions to be calculated
               const int end_index,
               const std::vector<int>& neighbour_list_start_index,
               const std::vector<int>& neighbour_list_end_index,
               const std::vector<int>& type_array, // type for atom
               const std::vector<int>& neighbour_list_array, // list of interactions between atom
               const std::vector<double>& x0_coord_array, // coord vectors for atoms
               const std::vector<double>& y0_coord_array,
               const std::vector<double>& z0_coord_array,
               std::vector<double>& x_coord_array, // coord vectors for atoms
               std::vector<double>& y_coord_array,
               std::vector<double>& z_coord_array,
               std::vector<double>& forces_array_x, //  vectors for forces
               std::vector<double>& forces_array_y,
               std::vector<double>& forces_array_z){

         return;
               }

     void stats_sld(){
      const int num_atoms=atoms::num_atoms;

      std::fill(sld::internal::fields_array_x.begin(), sld::internal::fields_array_x.end(), 0.0);
      std::fill(sld::internal::fields_array_y.begin(), sld::internal::fields_array_y.end(), 0.0);
      std::fill(sld::internal::fields_array_z.begin(), sld::internal::fields_array_z.end(), 0.0);
      std::fill(sld::internal::forces_array_x.begin(), sld::internal::forces_array_x.end(), 0.0);
      std::fill(sld::internal::forces_array_y.begin(), sld::internal::forces_array_y.end(), 0.0);
      std::fill(sld::internal::forces_array_z.begin(), sld::internal::forces_array_z.end(), 0.0);

      sld::compute_fields(0, // first atom for exchange interactions to be calculated
                        num_atoms, // last +1 atom to be calculated
                        atoms::neighbour_list_start_index,
                        atoms::neighbour_list_end_index,
                        atoms::type_array, // type for atom
                        atoms::neighbour_list_array, // list of interactions between atoms
                        atoms::x_coord_array,
                        atoms::y_coord_array,
                        atoms::z_coord_array,
                        atoms::x_spin_array,
                        atoms::y_spin_array,
                        atoms::z_spin_array,
                        sld::internal::forces_array_x,
                        sld::internal::forces_array_y,
                        sld::internal::forces_array_z,
                        sld::internal::fields_array_x,
                        sld::internal::fields_array_y,
                        sld::internal::fields_array_z);



      sld::J_eff= sld::compute_effective_J(0,atoms::num_atoms,
                  sld::internal::sumJ);
//
      sld::C_eff= sld::compute_effective_C(0,atoms::num_atoms,
                  sld::internal::sumC);



   }



   int suzuki_trotter(){
      const int num_atoms=atoms::num_atoms;
      sld::internal::prepare_sled_thermostat(); // prep the SLED thermostat for the current step
       // record coefficients after any SLED update so all ranks use the current thermostat bath temperature
      sld::internal::prepare_integrator_coefficients();
      double cay_dt=-mp::dt/4.0;//-dt4*consts::gyro - mp::dt contains gamma;
      double dt2=0.5*mp::dt_SI*1e12;



      // access the preallocated spin noise arrays rather than constructing three vectors on every rank and step
      std::vector<double>& Hx_th = sld::internal::spin_noise_array_x;
      std::vector<double>& Hy_th = sld::internal::spin_noise_array_y;
      std::vector<double>& Hz_th = sld::internal::spin_noise_array_z;
      // nonlinear field terms of the current single-spin visit, refilled by compute_spin_field
      static sld::internal::nonlinear_field_terms_t nonlinear_terms;

      generate (Hx_th.begin(),Hx_th.end(), mtrandom::gaussian);
      generate (Hy_th.begin(),Hy_th.end(), mtrandom::gaussian);
      generate (Hz_th.begin(),Hz_th.end(), mtrandom::gaussian);

      // access the equivalent preallocated lattice noise buffers
      std::vector<double>& Fx_th = sld::internal::lattice_noise_array_x;
      std::vector<double>& Fy_th = sld::internal::lattice_noise_array_y;
      std::vector<double>& Fz_th = sld::internal::lattice_noise_array_z;

      generate (Fx_th.begin(),Fx_th.end(), mtrandom::gaussian);
      generate (Fy_th.begin(),Fy_th.end(), mtrandom::gaussian);
      generate (Fz_th.begin(),Fz_th.end(), mtrandom::gaussian);


      std::fill(sld::internal::fields_array_x.begin(), sld::internal::fields_array_x.end(), 0.0);
      std::fill(sld::internal::fields_array_y.begin(), sld::internal::fields_array_y.end(), 0.0);
      std::fill(sld::internal::fields_array_z.begin(), sld::internal::fields_array_z.end(), 0.0);


      /*for(int atom=0;atom<=num_atoms-1;atom++){

         if ( (atom==555)) {
         std::cout<<std::setprecision(15)<<std::endl;

         std::cout <<" b i=atom="<<atom<<"\txyz "<<atoms::x_coord_array[atom]<<"\t"<<atoms::y_coord_array[atom]<<"\t"<<atoms::z_coord_array[atom]<<std::endl;
         std::cout <<" b i=atom="<<atom<<"\tvel "<<atoms::x_velo_array[atom]<<"\t"<<atoms::y_velo_array[atom]<<"\t"<<atoms::z_velo_array[atom]<<std::endl;


         }}*/

      for(int atom=0;atom<=num_atoms-1;atom++){




      sld::internal::compute_spin_field(atom, nonlinear_terms); // field and nonlinear terms only, for the single-spin update below

      sld::internal::update_single_spin(atom, cay_dt, Hx_th, Hy_th, Hz_th, nonlinear_terms);

      }


      std::fill(sld::internal::fields_array_x.begin(), sld::internal::fields_array_x.end(), 0.0);
      std::fill(sld::internal::fields_array_y.begin(), sld::internal::fields_array_y.end(), 0.0);
      std::fill(sld::internal::fields_array_z.begin(), sld::internal::fields_array_z.end(), 0.0);

      for(int atom=num_atoms-1;atom>=0;atom--){


         sld::internal::compute_spin_field(atom, nonlinear_terms); // field and nonlinear terms only, for the single-spin update below


         sld::internal::update_single_spin(atom, cay_dt, Hx_th, Hy_th, Hz_th, nonlinear_terms);

      }




      //forces are set to 0 for computation
      std::fill(sld::internal::forces_array_x.begin(), sld::internal::forces_array_x.end(), 0.0);
      std::fill(sld::internal::forces_array_y.begin(), sld::internal::forces_array_y.end(), 0.0);
      std::fill(sld::internal::forces_array_z.begin(), sld::internal::forces_array_z.end(), 0.0);



      sld::compute_fields(0, // first atom for exchange interactions to be calculated
                        num_atoms, // last +1 atom to be calculated
                        atoms::neighbour_list_start_index,
                        atoms::neighbour_list_end_index,
                        atoms::type_array, // type for atom
                        atoms::neighbour_list_array, // list of interactions between atoms
                        atoms::x_coord_array,
                        atoms::y_coord_array,
                        atoms::z_coord_array,
                        atoms::x_spin_array,
                        atoms::y_spin_array,
                        atoms::z_spin_array,
                        sld::internal::forces_array_x,
                        sld::internal::forces_array_y,
                        sld::internal::forces_array_z,
                        sld::internal::fields_array_x,
                        sld::internal::fields_array_y,
                        sld::internal::fields_array_z);

      sld::compute_forces(0, // first atom for exchange interactions to be calculated
                        num_atoms, // last +1 atom to be calculated
                        atoms::neighbour_list_start_index,
                        atoms::neighbour_list_end_index,
                        atoms::type_array, // type for atom
                        atoms::neighbour_list_array, // list of interactions between atoms
                        sld::internal::x0_coord_array, // list of isotropic exchange constants
                        sld::internal::y0_coord_array, // list of vectorial exchange constants
                        sld::internal::z0_coord_array, // list of tensorial exchange constants
                        atoms::x_coord_array,
                        atoms::y_coord_array,
                        atoms::z_coord_array,
                        sld::internal::forces_array_x,
                        sld::internal::forces_array_y,
                        sld::internal::forces_array_z,
                        sld::internal::potential_eng);


      internal::apply_phonon_pump(double(sim::time));

     //update position, Velocity
      for(int atom=0;atom<num_atoms;atom++){

      const unsigned int imat = atoms::type_array[atom];
      // look up the coefficients stored once for this material instead of recomputing them per atom
      const double dt2_m=sld::internal::lattice_dt2_over_mass_array[imat];
      const double f_eta=sld::internal::lattice_damping_factor_array[imat];
      const double velo_noise=sld::internal::lattice_noise_scale_array[imat];

             atoms::x_velo_array[atom] =  f_eta*atoms::x_velo_array[atom]+ dt2_m * sld::internal::forces_array_x[atom]+dt2*velo_noise*Fx_th[atom];
             atoms::y_velo_array[atom] =  f_eta*atoms::y_velo_array[atom]+ dt2_m * sld::internal::forces_array_y[atom]+dt2*velo_noise*Fy_th[atom];
             atoms::z_velo_array[atom] =  f_eta*atoms::z_velo_array[atom]+ dt2_m * sld::internal::forces_array_z[atom]+dt2*velo_noise*Fz_th[atom];

             atoms::x_coord_array[atom] +=  mp::dt_SI*1e12 * atoms::x_velo_array[atom];
             atoms::y_coord_array[atom] +=  mp::dt_SI*1e12 * atoms::y_velo_array[atom];
             atoms::z_coord_array[atom] +=  mp::dt_SI*1e12 * atoms::z_velo_array[atom];

       }




       internal::finish_phonon_pump_kick();

       // check positions vs max displacement
       sld::internal::check_displacement();

       //reset forces to 0 for v integration
        std::fill(sld::internal::forces_array_x.begin(), sld::internal::forces_array_x.end(), 0.0);
        std::fill(sld::internal::forces_array_y.begin(), sld::internal::forces_array_y.end(), 0.0);
        std::fill(sld::internal::forces_array_z.begin(), sld::internal::forces_array_z.end(), 0.0);


       sld::compute_fields(0, // first atom for exchange interactions to be calculated
                          num_atoms, // last +1 atom to be calculated
                          atoms::neighbour_list_start_index,
                          atoms::neighbour_list_end_index,
                          atoms::type_array, // type for atom
                          atoms::neighbour_list_array, // list of interactions between atoms
                          atoms::x_coord_array,
                          atoms::y_coord_array,
                          atoms::z_coord_array,
                          atoms::x_spin_array,
                          atoms::y_spin_array,
                          atoms::z_spin_array,
                          sld::internal::forces_array_x,
                          sld::internal::forces_array_y,
                          sld::internal::forces_array_z,
                          sld::internal::fields_array_x,
                          sld::internal::fields_array_y,
                          sld::internal::fields_array_z);


        sld::compute_forces(0, // first atom for exchange interactions to be calculated
                          num_atoms, // last +1 atom to be calculated
                          atoms::neighbour_list_start_index,
                          atoms::neighbour_list_end_index,
                          atoms::type_array, // type for atom
                          atoms::neighbour_list_array, // list of interactions between atoms
                          sld::internal::x0_coord_array, // list of isotropic exchange constants
                          sld::internal::y0_coord_array, // list of vectorial exchange constants
                          sld::internal::z0_coord_array, // list of tensorial exchange constants
                          atoms::x_coord_array,
                          atoms::y_coord_array,
                          atoms::z_coord_array,
                          sld::internal::forces_array_x,
                          sld::internal::forces_array_y,
                          sld::internal::forces_array_z,
                          sld::internal::potential_eng);


      internal::apply_phonon_pump(double(sim::time)+1.0);

      for(int atom=0;atom<num_atoms;atom++){

        const unsigned int imat = atoms::type_array[atom];
        // reuse the stored coefficients for the second velocity half step
        const double dt2_m=sld::internal::lattice_dt2_over_mass_array[imat];
        const double f_eta=sld::internal::lattice_damping_factor_array[imat];
        const double velo_noise=sld::internal::lattice_noise_scale_array[imat];


         atoms::x_velo_array[atom] =  f_eta*atoms::x_velo_array[atom] + dt2_m * sld::internal::forces_array_x[atom]+dt2*velo_noise*Fx_th[atom];
         atoms::y_velo_array[atom] =  f_eta*atoms::y_velo_array[atom] + dt2_m * sld::internal::forces_array_y[atom]+dt2*velo_noise*Fy_th[atom];
         atoms::z_velo_array[atom] =  f_eta*atoms::z_velo_array[atom] + dt2_m * sld::internal::forces_array_z[atom]+dt2*velo_noise*Fz_th[atom];


      }


      internal::finish_phonon_pump_kick();

  std::fill(sld::internal::fields_array_x.begin(), sld::internal::fields_array_x.end(), 0.0);
  std::fill(sld::internal::fields_array_y.begin(), sld::internal::fields_array_y.end(), 0.0);
  std::fill(sld::internal::fields_array_z.begin(), sld::internal::fields_array_z.end(), 0.0);


   for(int atom=0;atom<=num_atoms-1;atom++){


   sld::internal::compute_spin_field(atom, nonlinear_terms); // field and nonlinear terms only, for the single-spin update below

   sld::internal::update_single_spin(atom, cay_dt, Hx_th, Hy_th, Hz_th, nonlinear_terms);


   }

   std::fill(sld::internal::fields_array_x.begin(), sld::internal::fields_array_x.end(), 0.0);
   std::fill(sld::internal::fields_array_y.begin(), sld::internal::fields_array_y.end(), 0.0);
   std::fill(sld::internal::fields_array_z.begin(), sld::internal::fields_array_z.end(), 0.0);

   for(int atom=num_atoms-1;atom>=0;atom--){

      sld::internal::compute_spin_field(atom, nonlinear_terms); // field and nonlinear terms only, for the single-spin update below

      sld::internal::update_single_spin(atom, cay_dt, Hx_th, Hy_th, Hz_th, nonlinear_terms);

   }

    /*

      std::fill(sld::internal::fields_array_x.begin(), sld::internal::fields_array_x.end(), 0.0);
      std::fill(sld::internal::fields_array_y.begin(), sld::internal::fields_array_y.end(), 0.0);
      std::fill(sld::internal::fields_array_z.begin(), sld::internal::fields_array_z.end(), 0.0);

      sld::compute_fields(0, // first atom for exchange interactions to be calculated
                        num_atoms, // last +1 atom to be calculated
                        atoms::neighbour_list_start_index,
                        atoms::neighbour_list_end_index,
                        atoms::type_array, // type for atom
                        atoms::neighbour_list_array, // list of interactions between atoms
                        atoms::x_coord_array,
                        atoms::y_coord_array,
                        atoms::z_coord_array,
                        atoms::x_spin_array,
                        atoms::y_spin_array,
                        atoms::z_spin_array,
                        sld::internal::forces_array_x,
                        sld::internal::forces_array_y,
                        sld::internal::forces_array_z,
                        sld::internal::fields_array_x,
                        sld::internal::fields_array_y,
                        sld::internal::fields_array_z);

      sld::spin_temperature= sld::compute_spin_temperature(0,atoms::num_atoms,
                  atoms::type_array, // type for atom
                  atoms::x_spin_array,
                  atoms::y_spin_array,
                  atoms::z_spin_array,
                  sld::internal::fields_array_x,
                  sld::internal::fields_array_y,
                  sld::internal::fields_array_z,
                  mp::mu_s_array);
      //
      sld::lattice_temperature= sld::compute_lattice_temperature(0,atoms::num_atoms,
                  atoms::type_array,
                  atoms::x_velo_array,
                  atoms::y_velo_array,
                  atoms::z_velo_array);*/


               /*   for(int atom=0;atom<=num_atoms-1;atom++){

                     if ( (atom==555)) {
                     std::cout<<std::setprecision(15)<<std::endl;

                     std::cout <<" end i=atom="<<atom<<"\txyz "<<atoms::x_coord_array[atom]<<"\t"<<atoms::y_coord_array[atom]<<"\t"<<atoms::z_coord_array[atom]<<std::endl;
                     std::cout <<" end i=atom="<<atom<<"\tvel "<<atoms::x_velo_array[atom]<<"\t"<<atoms::y_velo_array[atom]<<"\t"<<atoms::z_velo_array[atom]<<std::endl;


                     }}*/

      sld::internal::update_sled_thermostat(); // update the SLED thermostat for the current step
      return EXIT_SUCCESS;
  }

namespace internal{

// Advance one spin over one single-spin update of duration dt/4.
// The field on spin i is B_i(s) = B_rest + B^nl_i(s), where B^nl_i includes all terms
// that depend on s_i itself (currently full Neel q, biquadratic exchange, second order
// uniaxial and fourth order cubic anisotropy). While positions and neighbour spins are fixed, B_rest is independent of s_i.
// B^nl_i is split into fixed-axis terms B_a, each constant during its own rotation, and composed symmetrically:
//    s_end = A_1(h/2) ... A_M(h/2) X(h) A_M(h/2) ... A_1(h/2) s_start
// where A_a rotates by eta B_a with eta = 1/(1+lambda^2) and X rotates by the central field C assembled below

void update_single_spin(const int atom,
                        const double cayley_dt,
                        const std::vector<double>& Hx_th,
                        const std::vector<double>& Hy_th,
                        const std::vector<double>& Hz_th,
                        const nonlinear_field_terms_t& nonlinear_terms){

   // When every term is linear in the spin being updated, don't carry out the extra decomposition since the field is fixed during the update
   if(!sld::internal::nonlinear_spin_hamiltonian){
      sld::internal::add_spin_noise(atom, atom+1, mp::dt_SI*1e12, atoms::type_array,
                                    atoms::x_spin_array, atoms::y_spin_array, atoms::z_spin_array,
                                    sld::internal::fields_array_x, sld::internal::fields_array_y, sld::internal::fields_array_z,
                                    Hx_th, Hy_th, Hz_th);
      sld::internal::cayley_update(atom, atom+1, cayley_dt,
                                   atoms::x_spin_array, atoms::y_spin_array, atoms::z_spin_array,
                                   sld::internal::fields_array_x, sld::internal::fields_array_y, sld::internal::fields_array_z);
      return;
   }

   const neel_vector_t spin_start = {atoms::x_spin_array[atom], atoms::y_spin_array[atom], atoms::z_spin_array[atom]};

   // The nonlinear field terms were collected by compute_spin_field for this visit, and positions and neighbour spins are fixed during it
   const neel_vector_t nonlinear_field_start = sld::internal::compute_nonlinear_field(nonlinear_terms, spin_start);
   static std::vector<spin_rotation_term_t> rotation_terms;
   sld::internal::build_spin_rotation_terms(nonlinear_terms, rotation_terms);

   // Central field C = eta [B_rest + N + D_0], frozen at the start of this single-spin visit (s_0 = spin_start):
   //    B_0 = B_rest + B^nl(s_0)   conservative field from compute_fields
   //    N = sigma xi               thermal field, drawn once per timestep and reused by every visit
   //    D_0 = lambda (B_0 + N) x s_0  damping correction of the complete field, so the nonlinear terms are still damped
   // The nonlinear precession eta B^nl(s) is applied by the rotations, so C + eta B^nl(s_0) is the field of the linear update (add_spin_noise)
   const unsigned int material = atoms::type_array[atom];
   const double lambda = sld::internal::spin_damping_array[material];
   const double eta = 1.0/(1.0+lambda*lambda);
   const double noise_scale = sld::internal::spin_noise_scale_array[material];
   const neel_vector_t conservative_field = {sld::internal::fields_array_x[atom], sld::internal::fields_array_y[atom], sld::internal::fields_array_z[atom]};
   const neel_vector_t remaining_field = neel_add(conservative_field, neel_scale(nonlinear_field_start, -1.0));
   const neel_vector_t thermal_field = {noise_scale*Hx_th[atom], noise_scale*Hy_th[atom], noise_scale*Hz_th[atom]};
   const neel_vector_t total_field = neel_add(conservative_field, thermal_field);
   const neel_vector_t total_cross_spin = {total_field.y*spin_start.z-total_field.z*spin_start.y, total_field.z*spin_start.x-total_field.x*spin_start.z, total_field.x*spin_start.y-total_field.y*spin_start.x};
   const neel_vector_t damping_correction = neel_scale(total_cross_spin, lambda);
   const neel_vector_t central_field = neel_scale(neel_add(neel_add(remaining_field, thermal_field), damping_correction), eta);

   // Half rotations about each fixed-axis term, a full rotation by the frozen central field, then the half rotations in reverse.
   const double half_dt = 0.5*cayley_dt*eta;
   neel_vector_t spin = spin_start;
   for(size_t n=0; n<rotation_terms.size(); ++n) spin = rotate_about_term_axis(rotation_terms[n], spin, half_dt);
   spin = sld::internal::cayley_update(spin, central_field, cayley_dt);
   for(size_t n=rotation_terms.size(); n-- > 0;) spin = rotate_about_term_axis(rotation_terms[n], spin, half_dt);

   atoms::x_spin_array[atom] = spin.x;
   atoms::y_spin_array[atom] = spin.y;
   atoms::z_spin_array[atom] = spin.z;
   return;
}

void cayley_update(const int start_index,
            const int end_index,
            double dt,
            std::vector<double>& x_spin_array, // coord vectors for atoms
            std::vector<double>& y_spin_array,
            std::vector<double>& z_spin_array,
            std::vector<double>& fields_array_x, //  vectors for fields
            std::vector<double>& fields_array_y,
            std::vector<double>& fields_array_z){

      for( int i = start_index; i<end_index; i++)
      {
          double Sx = x_spin_array[i];
          double Sy = y_spin_array[i];
          double Sz = z_spin_array[i];

          double Ax = fields_array_x[i] * dt;
          double Ay = fields_array_y[i] * dt;
          double Az = fields_array_z[i] * dt;

          double AS = Ax*Sx + Ay*Sy + Az*Sz;
          double A2 = Ax * Ax + Ay* Ay + Az * Az;

          double AxSx = Ay * Sz - Az * Sy;
          double AxSy = Az * Sx - Ax * Sz;
          double AxSz = Ax * Sy - Ay * Sx;

          double factor = 1.0 / (1.0 + 0.25 * A2);

          x_spin_array[i] = (Sx * ( 1.0 - 0.25 * A2) + AxSx + 0.5 * Ax * AS) * factor;
          y_spin_array[i] = (Sy * ( 1.0 - 0.25 * A2) + AxSy + 0.5 * Ay * AS) * factor;
          z_spin_array[i] = (Sz * ( 1.0 - 0.25 * A2) + AxSz + 0.5 * Az * AS) * factor;
      }

  return;
}

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
            const std::vector<double>& Hx_th, // vectors for fields
            const std::vector<double>& Hy_th,
            const std::vector<double>& Hz_th){


     for( int i = start_index; i<end_index; i++)

    {
        // use the material index passed with this spin array
        const unsigned int imat = type_array[i];
        // these values are already filtered based on equilibration or production stages
        const double lambda=sld::internal::spin_damping_array[imat];
        const double spin_noise=sld::internal::spin_noise_scale_array[imat];



        double Sx = x_spin_array[i];
        double Sy = y_spin_array[i];
        double Sz = z_spin_array[i];

        double Fx = fields_array_x[i] + spin_noise * Hx_th[i];
        double Fy = fields_array_y[i] + spin_noise * Hy_th[i];
        double Fz = fields_array_z[i] + spin_noise * Hz_th[i];

        double FxSx = Fy * Sz - Fz * Sy;
        double FxSy = Fz * Sx - Fx * Sz;
        double FxSz = Fx * Sy - Fy * Sx;

        double inv_l2 = 1.0 / (1.0 + lambda*lambda);

        fields_array_x[i] = (Fx + lambda * FxSx) * inv_l2;
        fields_array_y[i] = (Fy + lambda * FxSy) * inv_l2;
        fields_array_z[i] = (Fz + lambda * FxSz) * inv_l2;
    }

return;
}//end of add_spin_noise

// Rotate one trial spin by the Cayley map, the same map as above, returning the result instead of writing it to the atom
neel_vector_t cayley_update(const neel_vector_t& spin,
                            const neel_vector_t& field,
                            const double dt){

      const neel_vector_t A = neel_scale(field, dt);
      const double AS = neel_dot(A, spin);
      const double A2 = neel_dot(A, A);
      const neel_vector_t cross = {A.y*spin.z-A.z*spin.y, A.z*spin.x-A.x*spin.z, A.x*spin.y-A.y*spin.x};
      const double factor = 1.0/(1.0+0.25*A2);
      return neel_vector_t{(spin.x*(1.0-0.25*A2)+cross.x+0.5*A.x*AS)*factor,
                           (spin.y*(1.0-0.25*A2)+cross.y+0.5*A.y*AS)*factor,
                           (spin.z*(1.0-0.25*A2)+cross.z+0.5*A.z*AS)*factor};
}


} // end of internal namespace

} // end of sld namespace
