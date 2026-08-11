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
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

// Vampire headers
#include "atoms.hpp"
#include "constants.hpp"
#include "create.hpp"
#include "material.hpp"
#include "sld.hpp"
#include "vmpi.hpp"

// sld module headers
#include "internal.hpp"


namespace sld{

namespace{

// solve I omega = L for the inertia tensor {Ixx, Iyy, Izz, Ixy, Ixz, Iyz}
bool calculate_angular_velocity(const double inertia[6],
                                const double angular_momentum[3],
                                double angular_velocity[3]){

   // I components
   const double ixx = inertia[0];
   const double iyy = inertia[1];
   const double izz = inertia[2];
   const double ixy = inertia[3];
   const double ixz = inertia[4];
   const double iyz = inertia[5];

   const double determinant = ixx*(iyy*izz-iyz*iyz) - ixy*(ixy*izz-ixz*iyz) + ixz*(ixy*iyz-ixz*iyy);
   // scale by the largest component of tensor to avoid too small determinant
   const double scale = std::max(std::max(std::fabs(ixx), std::fabs(iyy)),
                                 std::max(std::fabs(izz),
                                 std::max(std::fabs(ixy),
                                 std::max(std::fabs(ixz), std::fabs(iyz)))));
   if(scale <= 0.0 || std::fabs(determinant) <= 64.0*std::numeric_limits<double>::epsilon()*scale*scale*scale){
      return false;
   }

   const double inverse_determinant = 1.0/determinant;
   const double inverse_xx = (iyy*izz-iyz*iyz)*inverse_determinant;
   const double inverse_xy = (ixz*iyz-ixy*izz)*inverse_determinant;
   const double inverse_xz = (ixy*iyz-ixz*iyy)*inverse_determinant;
   const double inverse_yy = (ixx*izz-ixz*ixz)*inverse_determinant;
   const double inverse_yz = (ixy*ixz-ixx*iyz)*inverse_determinant;
   const double inverse_zz = (ixx*iyy-ixy*ixy)*inverse_determinant;

   angular_velocity[0] = inverse_xx*angular_momentum[0] +
                         inverse_xy*angular_momentum[1] +
                         inverse_xz*angular_momentum[2];
   angular_velocity[1] = inverse_xy*angular_momentum[0] +
                         inverse_yy*angular_momentum[1] +
                         inverse_yz*angular_momentum[2];
   angular_velocity[2] = inverse_xz*angular_momentum[0] +
                         inverse_yz*angular_momentum[1] +
                         inverse_zz*angular_momentum[2];
   return true;
}

} // end of anonymous namespace


   // Calculate the spin temperature from fields
   double compute_spin_temperature(const int start_index, // first atom for exchange interactions to be calculated
               const int end_index,
               const std::vector<int>& type_array, // type for atom
               const std::vector<double>& x_spin_array, // coord vectors for atoms
               const std::vector<double>& y_spin_array,
               const std::vector<double>& z_spin_array,
               const std::vector<double>& fields_array_x, //  vectors for fields
               const std::vector<double>& fields_array_y,
               const std::vector<double>& fields_array_z,
               const std::vector <double>& mu_s_array){

                double SxH2=0.0;
                double denominator=0.0;
                for (int at=start_index;at<end_index;at++){
                    double Sx = x_spin_array[at];
                    double Sy = y_spin_array[at];
                    double Sz = z_spin_array[at];
                    double Hx = fields_array_x[at];
                    double Hy = fields_array_y[at];
                    double Hz = fields_array_z[at];

                     double SxHx = Sy * Hz - Sz * Hy;
                     double SxHy = Sz * Hx - Sx * Hz;
                     double SxHz = Sx * Hy - Sy * Hx;
                     // weight each atom by its own material moment
                     SxH2 += mu_s_array[type_array[at]] * (SxHx*SxHx + SxHy*SxHy + SxHz*SxHz);
                     // D_i=2 s_i.H_i-Tr[(I-s_i s_i^T)dH_i/ds_i].  The trace is zero for Hamiltonians linear in each spin and is accumulated for supported nonlinear terms including full neel coupling, second order uniaxial anisotropy and fourth order cubic anisotropy
                     denominator += 2.0*(Sx*Hx + Sy*Hy + Sz*Hz);
                     if(sld::internal::spin_temperature_correction){
                        denominator -= sld::internal::spin_hessian_trace[at]; // deduct the nonlinear term by defualt
                     }

                }
               #ifdef MPICF
                  double sums[2] = {SxH2, denominator};
                  MPI_Allreduce(MPI_IN_PLACE, sums, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                  SxH2 = sums[0];
                  denominator = sums[1];
               #endif

               if(denominator <= 0.0) return 0.0; 
               // T_s=mu_B Sum_i m_i|s_i x H_i|^2/(k_B Sum_i D_i)
               double T_spin = constants::muB / constants::kB * SxH2 / denominator;

      return T_spin;

      }//end of spin_temperature


double compute_lattice_temperature(const int start_index, // first atom for exchange interactions to be calculated
            const int end_index,
            const std::vector<int>& type_array, // type for atom
            const std::vector<double>& velo_array_x, // coord vectors for atoms
            const std::vector<double>& velo_array_y,
            const std::vector<double>& velo_array_z){

   double mass = 0.0;
   double momentum_x = 0.0;
   double momentum_y = 0.0;
   double momentum_z = 0.0;
   double mass_position_x = 0.0;
   double mass_position_y = 0.0;
   double mass_position_z = 0.0;
   double mass_velocity_squared = 0.0;
   double number_of_atoms = 0.0;

   // Sum M, P = Sum_i m_i v_i, M r_cm = Sum_i m_i r_i and Sum_i m_i |v_i|^2 and m_i is the material specific mass
   for(int at=start_index; at<end_index; ++at){
      const double atom_mass = sld::internal::mp[type_array[at]].mass.get();
      const double vx = velo_array_x[at];
      const double vy = velo_array_y[at];
      const double vz = velo_array_z[at];

      mass += atom_mass;
      momentum_x += atom_mass * vx;
      momentum_y += atom_mass * vy;
      momentum_z += atom_mass * vz;
      mass_position_x += atom_mass * atoms::x_coord_array[at];
      mass_position_y += atom_mass * atoms::y_coord_array[at];
      mass_position_z += atom_mass * atoms::z_coord_array[at];
      mass_velocity_squared += atom_mass * (vx*vx + vy*vy + vz*vz);
      number_of_atoms += 1.0;
   }

   // construct global sums before applying either rigid-body correction
   #ifdef MPICF
      double sums[9] = {mass, momentum_x, momentum_y, momentum_z,
                        mass_position_x, mass_position_y, mass_position_z,
                        mass_velocity_squared, number_of_atoms};
      MPI_Allreduce(MPI_IN_PLACE, sums, 9, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      mass = sums[0];
      momentum_x = sums[1];
      momentum_y = sums[2];
      momentum_z = sums[3];
      mass_position_x = sums[4];
      mass_position_y = sums[5];
      mass_position_z = sums[6];
      mass_velocity_squared = sums[7];
      number_of_atoms = sums[8];
   #endif

   if(mass <= 0.0 || number_of_atoms <= 0.0) return 0.0;

   // Sum_i m_i |v_i-v_cm|^2 = Sum_i m_i |v_i|^2-|Sum_i m_i v_i|^2/Sum_i m_i
   const double centre_of_mass_term = (momentum_x*momentum_x + momentum_y*momentum_y + momentum_z*momentum_z) / mass;
   double thermal_mass_velocity_squared = std::max(0.0, mass_velocity_squared - centre_of_mass_term);

   // removing COM translation removes 3 DOFs
   double degrees_of_freedom = 3.0*number_of_atoms-3.0;

   // if any direction is periodic, rotational subtraction is disabled
   const bool remove_rigid_body_rotation = !cs::pbc[0] && !cs::pbc[1] && !cs::pbc[2];
   if(remove_rigid_body_rotation && degrees_of_freedom >= 3.0){
      const double centre_of_mass_x = mass_position_x/mass;
      const double centre_of_mass_y = mass_position_y/mass;
      const double centre_of_mass_z = mass_position_z/mass;
      const double centre_of_mass_velocity_x = momentum_x/mass;
      const double centre_of_mass_velocity_y = momentum_y/mass;
      const double centre_of_mass_velocity_z = momentum_z/mass;

      double inertia[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
      double angular_momentum[3] = {0.0, 0.0, 0.0};
      for(int at=start_index; at<end_index; ++at){
         const double atom_mass = sld::internal::mp[type_array[at]].mass.get();
         const double rx = atoms::x_coord_array[at]-centre_of_mass_x;
         const double ry = atoms::y_coord_array[at]-centre_of_mass_y;
         const double rz = atoms::z_coord_array[at]-centre_of_mass_z;
         const double vx = velo_array_x[at]-centre_of_mass_velocity_x;
         const double vy = velo_array_y[at]-centre_of_mass_velocity_y;
         const double vz = velo_array_z[at]-centre_of_mass_velocity_z;

         // I = Sum_i m_i [(r_i.r_i) 1-r_i r_i^T]
         inertia[0] += atom_mass*(ry*ry+rz*rz); // Ixx
         inertia[1] += atom_mass*(rx*rx+rz*rz); // Iyy
         inertia[2] += atom_mass*(rx*rx+ry*ry); // Izz
         inertia[3] -= atom_mass*rx*ry;          // Ixy
         inertia[4] -= atom_mass*rx*rz;          // Ixz
         inertia[5] -= atom_mass*ry*rz;          // Iyz

         // L = Sum_i (r_i-r_cm) x m_i(v_i-v_cm)
         angular_momentum[0] += atom_mass*(ry*vz-rz*vy);
         angular_momentum[1] += atom_mass*(rz*vx-rx*vz);
         angular_momentum[2] += atom_mass*(rx*vy-ry*vx);
      }

      #ifdef MPICF
         MPI_Allreduce(MPI_IN_PLACE, inertia, 6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, angular_momentum, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif

      double angular_velocity[3] = {0.0, 0.0, 0.0};
      if(calculate_angular_velocity(inertia, angular_momentum, angular_velocity)){
         // Sum_i m_i |v_rot,i|^2 = omega.L, where v_rot,i = omega x (r_i-r_cm) and omega = I^-1 L
         const double rotational_mass_velocity_squared = angular_velocity[0]*angular_momentum[0] + angular_velocity[1]*angular_momentum[1] + angular_velocity[2]*angular_momentum[2];
         thermal_mass_velocity_squared = std::max(0.0, thermal_mass_velocity_squared-rotational_mass_velocity_squared);
         degrees_of_freedom -= 3.0;
      }
   }

   if(degrees_of_freedom <= 0.0) return 0.0;

   // T_l = Sum_i m_i |v_i-v_cm-v_rot,i|^2/(f k_B)
   const double T_lat = thermal_mass_velocity_squared / (degrees_of_freedom * constants::kB_eV);

   return T_lat;

}//end of lattice_temperature

   } // end of sld namespace
