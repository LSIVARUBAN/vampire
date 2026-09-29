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
#include <sstream>
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




      sld::compute_fields(atom, // first atom for exchange interactions to be calculated
                        atom+1, // last +1 atom to be calculated
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

      sld::internal::update_single_spin(atom, cay_dt, Hx_th, Hy_th, Hz_th);

      }


      std::fill(sld::internal::fields_array_x.begin(), sld::internal::fields_array_x.end(), 0.0);
      std::fill(sld::internal::fields_array_y.begin(), sld::internal::fields_array_y.end(), 0.0);
      std::fill(sld::internal::fields_array_z.begin(), sld::internal::fields_array_z.end(), 0.0);

      for(int atom=num_atoms-1;atom>=0;atom--){


         sld::compute_fields(atom, // first atom for exchange interactions to be calculated
                           atom+1, // last +1 atom to be calculated
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


         sld::internal::update_single_spin(atom, cay_dt, Hx_th, Hy_th, Hz_th);

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


  std::fill(sld::internal::fields_array_x.begin(), sld::internal::fields_array_x.end(), 0.0);
  std::fill(sld::internal::fields_array_y.begin(), sld::internal::fields_array_y.end(), 0.0);
  std::fill(sld::internal::fields_array_z.begin(), sld::internal::fields_array_z.end(), 0.0);


   for(int atom=0;atom<=num_atoms-1;atom++){


   sld::compute_fields(atom, // first atom for exchange interactions to be calculated
                     atom+1, // last +1 atom to be calculated
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

   sld::internal::update_single_spin(atom, cay_dt, Hx_th, Hy_th, Hz_th);


   }

   std::fill(sld::internal::fields_array_x.begin(), sld::internal::fields_array_x.end(), 0.0);
   std::fill(sld::internal::fields_array_y.begin(), sld::internal::fields_array_y.end(), 0.0);
   std::fill(sld::internal::fields_array_z.begin(), sld::internal::fields_array_z.end(), 0.0);

   for(int atom=num_atoms-1;atom>=0;atom--){

      sld::compute_fields(atom, // first atom for exchange interactions to be calculated
                        atom+1, // last +1 atom to be calculated
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

      sld::internal::update_single_spin(atom, cay_dt, Hx_th, Hy_th, Hz_th);

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

// Advance one spin, accepting the nonlinear midpoint when its iteration converges.
// The field on spin i is B_i(s) = B_rest + B^nl_i(s), where B^nl_i includes all terms
// that depend on s_i itself (currrently supports full Neel q, biquadratic exchange, second order
// uniaxial and fourth order cubic anisotropy). While positions and neighbour spins are fixed, B_rest is independent of s_i
// The update solves
//    s_end - s_start = h D(m, B_rest + B^nl_i(m)) x m,   m = (s_start + s_end)/2,
// D applies damping and noise and h = cayley_dt
void update_single_spin(const int atom,
                        const double cayley_dt,
                        const std::vector<double>& Hx_th,
                        const std::vector<double>& Hy_th,
                        const std::vector<double>& Hz_th){

   // When every term is linear in the spin being updated keep the original approach since the fixed field is unchanged during the update
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

   // Gather the nonlinear field terms once. Positions and neighbour spins are fixed during this update,
   // so each trial spin reuses the same nonlinear field terms instead of recalculating them.
   static nonlinear_field_terms_t nonlinear_terms;
   sld::internal::gather_nonlinear_field_terms(atom, nonlinear_terms);
   const neel_vector_t nonlinear_start = sld::internal::compute_nonlinear_field(nonlinear_terms, spin_start);

   // B_rest = B(s_start) - B^nl(s_start) includes every field term that is independent of s_i
   const neel_vector_t field_rest = {sld::internal::fields_array_x[atom]-nonlinear_start.x,
                                     sld::internal::fields_array_y[atom]-nonlinear_start.y,
                                     sld::internal::fields_array_z[atom]-nonlinear_start.z};

   // Use the original fixed-field result as an initial guess, then solve the Cayley midpoint equation by fixed-point iteration
   neel_vector_t effective = sld::internal::add_spin_noise(atom, spin_start, neel_add(field_rest, nonlinear_start), Hx_th, Hy_th, Hz_th);
   neel_vector_t candidate = sld::internal::cayley_update(spin_start, effective, cayley_dt);

   const int maximum_iterations = 50; // maximum number of iterations to attempt before aborting the simulation
   const double tolerance_squared = 1.0e-24; // squared tolerance for convergence of the midpoint spin update
   bool converged = false;
   double change_squared = 0.0;
   for(int iteration=0; iteration<maximum_iterations; ++iteration){
      const neel_vector_t midpoint = neel_scale(neel_add(spin_start, candidate), 0.5);
      const neel_vector_t nonlinear_midpoint = sld::internal::compute_nonlinear_field(nonlinear_terms, midpoint);
      effective = sld::internal::add_spin_noise(atom, midpoint, neel_add(field_rest, nonlinear_midpoint), Hx_th, Hy_th, Hz_th);
      const neel_vector_t next = sld::internal::cayley_update(spin_start, effective, cayley_dt);
      const neel_vector_t change = neel_add(next, neel_scale(candidate, -1.0));
      change_squared = neel_dot(change, change);
      if(!std::isfinite(change_squared)){
         std::ostringstream message;
         message << "Non-finite midpoint iteration for local atom " << atom << " at step " << sim::time << ". Aborting simulation.";
         err::zexit(message.str());
      }
      candidate = next;
      if(change_squared <= tolerance_squared){
         converged = true;
         break;
      }
   }

   if(!converged){
      std::ostringstream message;
      message << "Midpoint spin solve failed for atom " << atom << " at step " << sim::time << ". Aborting simulation.";
      err::zexit(message.str());
   }

   atoms::x_spin_array[atom] = candidate.x;
   atoms::y_spin_array[atom] = candidate.y;
   atoms::z_spin_array[atom] = candidate.z;
   sld::internal::fields_array_x[atom] = effective.x;
   sld::internal::fields_array_y[atom] = effective.y;
   sld::internal::fields_array_z[atom] = effective.z;
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

// Rotate one trial spin by the Cayley map using the same using the same math as above but return the result instead of writing it to the atom
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

// Apply the damping and spin noise as above just to one field using a supplied trial spin instead of the stored spin
neel_vector_t add_spin_noise(const int atom,
                             const neel_vector_t& spin,
                             const neel_vector_t& field,
                             const std::vector<double>& Hx_th,
                             const std::vector<double>& Hy_th,
                             const std::vector<double>& Hz_th){
   
      const unsigned int imat = atoms::type_array[atom];
      const double lambda = sld::internal::spin_damping_array[imat];
      const double noise = sld::internal::spin_noise_scale_array[imat];
      const neel_vector_t F = {field.x+noise*Hx_th[atom], field.y+noise*Hy_th[atom], field.z+noise*Hz_th[atom]};
      const neel_vector_t FxS = {F.y*spin.z-F.z*spin.y, F.z*spin.x-F.x*spin.z, F.x*spin.y-F.y*spin.x};
      return neel_scale(neel_add(F, neel_scale(FxS, lambda)), 1.0/(1.0+lambda*lambda));
}


} // end of internal namespace

} // end of sld namespace
