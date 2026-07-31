//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) Mara Strungaru 2022. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <limits>
#include <string>

// Vampire headers
#include "constants.hpp"
#include "create.hpp"
#include "errors.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vmpi.hpp"
#include "vio.hpp"
#include "sld.hpp"
#include "atoms.hpp"


namespace stats{

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


//------------------------------------------------------------------------------------------------------
// Function to determine if class is properly initialized
//------------------------------------------------------------------------------------------------------
bool lattice_temp_statistic_t::is_initialized(){
   return initialized;
}

//------------------------------------------------------------------------------------------------------
// Function to initialize mask
//------------------------------------------------------------------------------------------------------
void lattice_temp_statistic_t::set_mask(const int in_mask_size, std::vector<int> in_mask, const std::vector<double>& mm){

   // Check that mask values never exceed mask_size
   for(unsigned int atom=0; atom<in_mask.size(); ++atom){
      if(in_mask[atom] > in_mask_size-1){
         terminaltextcolor(RED);
         std::cerr << "Programmer Error - mask id " << in_mask[atom] << " is greater than number of elements for mask "<< in_mask_size << std::endl;
         terminaltextcolor(WHITE);
         zlog << zTs() << "Programmer Error - mask id " << in_mask[atom] << " is greater than number of elements for mask "<< in_mask_size << std::endl;
         err::vexit();
      }
   }

   // save mask to internal storage
   num_atoms = in_mask.size();
   mask_size = in_mask_size - 1; // last element contains energy for non-magnetic atoms
   mean_counter = 0.0;
   mask=in_mask; // copy contents of vector
   lattice_temp.resize(in_mask_size, 0.0);
   mean_lattice_temp.resize(in_mask_size, 0.0);
   normalisation.resize(in_mask_size, 0.0);




   // calculate number of spins in each mask
   for(int atom=0; atom<num_atoms; ++atom){
      const int mask_id = mask[atom]; // get mask id
      normalisation[mask_id] += 1.0;
   }

   // Calculate normalisation for all CPUs
   #ifdef MPICF
      MPI_Allreduce(MPI_IN_PLACE, &normalisation[0], mask_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   #endif



   // determine mask id's with no atoms
   num_atoms_in_mask.resize(in_mask_size,0);
   for(unsigned int atom=0; atom<in_mask.size(); ++atom){
      int mask_id = in_mask[atom];
      // add atoms to mask
      num_atoms_in_mask[mask_id]++;
   }

   // Reduce on all CPUs
   #ifdef MPICF
      MPI_Allreduce(MPI_IN_PLACE, &num_atoms_in_mask[0], mask_size, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   #endif

   // Check for no atoms in mask on any CPU
   for(int mask_id=0; mask_id<mask_size; ++mask_id){
      // if no atoms exist then add to zero list
      if(num_atoms_in_mask[mask_id]==0){
         zero_list.push_back(mask_id);
      }
   }

   // Set flag indicating correct initialization
   initialized = true;

   return;

}

//------------------------------------------------------------------------------------------------------
// Function to get mask needed for gpu acceleration of statistics calculation
//------------------------------------------------------------------------------------------------------
void lattice_temp_statistic_t::get_mask(std::vector<int>& out_mask,  std::vector<double>& out_normalisation){

   // copy data to objects
   out_mask = mask;
   out_normalisation = normalisation;


   return;

}


void lattice_temp_statistic_t::calculate_lattice_temp(const std::vector<double>& velo_array_x, // coord vectors for atoms
            const std::vector<double>& velo_array_y,
            const std::vector<double>& velo_array_z){

   std::fill(lattice_temp.begin(),lattice_temp.end(),0.0);
   std::vector<double> mass(mask_size, 0.0);
   std::vector<double> momentum_x(mask_size, 0.0);
   std::vector<double> momentum_y(mask_size, 0.0);
   std::vector<double> momentum_z(mask_size, 0.0);
   std::vector<double> mass_position_x(mask_size, 0.0);
   std::vector<double> mass_position_y(mask_size, 0.0);
   std::vector<double> mass_position_z(mask_size, 0.0);

   // Sum M, P = Sum_i m_i v_i, M r_cm = Sum_i m_i r_i and Sum_i m_i |v_i|^2 and m_i is the material specific mass
   for(int atom=0; atom < num_atoms; ++atom){

      const unsigned int mask_id = mask[atom]; // get mask id
      if(mask_id >= static_cast<unsigned int>(mask_size)) continue;

      const double atom_mass = atoms::mass_spin_array[atom];
      const double vx = velo_array_x[atom];
      const double vy = velo_array_y[atom];
      const double vz = velo_array_z[atom];

      mass[mask_id] += atom_mass;
      momentum_x[mask_id] += atom_mass * vx;
      momentum_y[mask_id] += atom_mass * vy;
      momentum_z[mask_id] += atom_mass * vz;
      mass_position_x[mask_id] += atom_mass * atoms::x_coord_array[atom];
      mass_position_y[mask_id] += atom_mass * atoms::y_coord_array[atom];
      mass_position_z[mask_id] += atom_mass * atoms::z_coord_array[atom];
      lattice_temp[mask_id] += atom_mass * (vx*vx + vy*vy + vz*vz);
   }

   // construct global sums before applying either rigid-body correction
   #ifdef MPICF
      MPI_Allreduce(MPI_IN_PLACE, &lattice_temp[0], mask_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &mass[0], mask_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &momentum_x[0], mask_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &momentum_y[0], mask_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &momentum_z[0], mask_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &mass_position_x[0], mask_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &mass_position_y[0], mask_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &mass_position_z[0], mask_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   #endif

   // if any direction is periodic, rotational subtraction is disabled
   const bool remove_rigid_body_rotation = !cs::pbc[0] && !cs::pbc[1] && !cs::pbc[2];
   std::vector<double> inertia_xx(mask_size, 0.0);
   std::vector<double> inertia_yy(mask_size, 0.0);
   std::vector<double> inertia_zz(mask_size, 0.0);
   std::vector<double> inertia_xy(mask_size, 0.0);
   std::vector<double> inertia_xz(mask_size, 0.0);
   std::vector<double> inertia_yz(mask_size, 0.0);
   std::vector<double> angular_momentum_x(mask_size, 0.0);
   std::vector<double> angular_momentum_y(mask_size, 0.0);
   std::vector<double> angular_momentum_z(mask_size, 0.0);

   if(remove_rigid_body_rotation){
      for(int atom=0; atom<num_atoms; ++atom){
         const unsigned int mask_id = mask[atom];
         if(mask_id >= static_cast<unsigned int>(mask_size) || mass[mask_id] <= 0.0) continue;

         const double atom_mass = atoms::mass_spin_array[atom];
         const double centre_of_mass_x = mass_position_x[mask_id]/mass[mask_id];
         const double centre_of_mass_y = mass_position_y[mask_id]/mass[mask_id];
         const double centre_of_mass_z = mass_position_z[mask_id]/mass[mask_id];
         const double centre_of_mass_velocity_x = momentum_x[mask_id]/mass[mask_id];
         const double centre_of_mass_velocity_y = momentum_y[mask_id]/mass[mask_id];
         const double centre_of_mass_velocity_z = momentum_z[mask_id]/mass[mask_id];
         const double rx = atoms::x_coord_array[atom]-centre_of_mass_x;
         const double ry = atoms::y_coord_array[atom]-centre_of_mass_y;
         const double rz = atoms::z_coord_array[atom]-centre_of_mass_z;
         const double vx = velo_array_x[atom]-centre_of_mass_velocity_x;
         const double vy = velo_array_y[atom]-centre_of_mass_velocity_y;
         const double vz = velo_array_z[atom]-centre_of_mass_velocity_z;

         // I = Sum_i m_i [(r_i.r_i) 1-r_i r_i^T]
         inertia_xx[mask_id] += atom_mass*(ry*ry+rz*rz);
         inertia_yy[mask_id] += atom_mass*(rx*rx+rz*rz);
         inertia_zz[mask_id] += atom_mass*(rx*rx+ry*ry);
         inertia_xy[mask_id] -= atom_mass*rx*ry;
         inertia_xz[mask_id] -= atom_mass*rx*rz;
         inertia_yz[mask_id] -= atom_mass*ry*rz;

         // L = Sum_i (r_i-r_cm) x m_i(v_i-v_cm)
         angular_momentum_x[mask_id] += atom_mass*(ry*vz-rz*vy);
         angular_momentum_y[mask_id] += atom_mass*(rz*vx-rx*vz);
         angular_momentum_z[mask_id] += atom_mass*(rx*vy-ry*vx);
      }

      #ifdef MPICF
         MPI_Allreduce(MPI_IN_PLACE, &inertia_xx[0], mask_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, &inertia_yy[0], mask_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, &inertia_zz[0], mask_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, &inertia_xy[0], mask_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, &inertia_xz[0], mask_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, &inertia_yz[0], mask_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, &angular_momentum_x[0], mask_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, &angular_momentum_y[0], mask_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, &angular_momentum_z[0], mask_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif
   }

   for(int mask_id=0; mask_id < mask_size; ++mask_id){
      if(num_atoms_in_mask[mask_id] <= 0 || mass[mask_id] <= 0.0){
         lattice_temp[mask_id] = 0.0;
         continue;
      }

      // Sum_i m_i |v_i-v_cm|^2 = Sum_i m_i |v_i|^2-|Sum_i m_i v_i|^2/Sum_i m_i
      const double centre_of_mass_term = (momentum_x[mask_id]*momentum_x[mask_id] + momentum_y[mask_id]*momentum_y[mask_id] + momentum_z[mask_id]*momentum_z[mask_id]) / mass[mask_id];
      double corrected_mass_velocity_squared = std::max(0.0, lattice_temp[mask_id] - centre_of_mass_term);

      // removing COM translation removes 3 DOFs
      double degrees_of_freedom = 3.0*double(num_atoms_in_mask[mask_id])-3.0;

      if(remove_rigid_body_rotation && degrees_of_freedom >= 3.0){
         const double inertia[6] = {inertia_xx[mask_id], inertia_yy[mask_id],
                                    inertia_zz[mask_id], inertia_xy[mask_id],
                                    inertia_xz[mask_id], inertia_yz[mask_id]};
         const double angular_momentum[3] = {angular_momentum_x[mask_id],
                                             angular_momentum_y[mask_id],
                                             angular_momentum_z[mask_id]};
         double angular_velocity[3] = {0.0, 0.0, 0.0};
         if(calculate_angular_velocity(inertia, angular_momentum, angular_velocity)){
            // Sum_i m_i |v_rot,i|^2 = omega.L, where v_rot,i = omega x (r_i-r_cm) and omega = I^-1 L
            const double rotational_mass_velocity_squared = angular_velocity[0]*angular_momentum[0] + angular_velocity[1]*angular_momentum[1] + angular_velocity[2]*angular_momentum[2];
            corrected_mass_velocity_squared = std::max(0.0, corrected_mass_velocity_squared-rotational_mass_velocity_squared);
            degrees_of_freedom -= 3.0;
         }
      }

      if(degrees_of_freedom <= 0.0){
         lattice_temp[mask_id] = 0.0;
         continue;
      }

      // store 3 Sum_i m_i |v_i-v_cm-v_rot,i|^2/f so the existing output conversion by 3 k_B gives T_l
      lattice_temp[mask_id] = 3.0*corrected_mass_velocity_squared/degrees_of_freedom;
   }


   // Zero empty mask id's
   for(unsigned int id=0; id<zero_list.size(); ++id) lattice_temp[zero_list[id]]=0.0;

   const int tsize = lattice_temp.size();
   for(int idx = 0; idx < tsize; ++idx) mean_lattice_temp[idx] += lattice_temp[idx];
   mean_counter+=1.0;
   // std::cout<<"lat2"<<std::endl;


   return;

}

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
const std::vector<double>& lattice_temp_statistic_t::get_lattice_temp(){

   return lattice_temp;

}

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
void lattice_temp_statistic_t::set_lattice_temp(std::vector<double>& new_lattice_temp, std::vector<double>& new_mean_lattice_temp, long counter){

   lattice_temp = new_lattice_temp;

   const size_t array_size = mean_lattice_temp.size();
   for(size_t i=0; i< array_size; ++i){
      mean_lattice_temp[i] += new_mean_lattice_temp[i];
   }

   // update counter
   mean_counter += double(counter);

}

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
void lattice_temp_statistic_t::reset_lattice_temp_averages(){

   std::fill(mean_lattice_temp.begin(),mean_lattice_temp.end(),0.0);

   // reset data counter
   mean_counter = 0.0;

   return;

}

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
std::string lattice_temp_statistic_t::output_lattice_temp(bool header){

   // result string stream
   std::ostringstream res;

   // set custom precision if enabled
   if(vout::custom_precision){
      res.precision(vout::precision);
      if(vout::fixed) res.setf( std::ios::fixed, std::ios::floatfield );
   }
   vout::fixed_width_output result(res,vout::fw_size);

   for(int mask_id=0; mask_id<mask_size; ++mask_id){
      if(header){
         result << name + std::to_string(mask_id) + "_TL";
      }
      else{
         result << 0.5*( 2.0 /(3.0*constants::kB_eV)) * lattice_temp[mask_id ];
      }
   }

   return result.str();

}

//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
std::string lattice_temp_statistic_t::output_mean_lattice_temp(bool header){

   // result string stream
   std::ostringstream res;

   // set custom precision if enabled
   if(vout::custom_precision){
      res.precision(vout::precision);
      if(vout::fixed) res.setf( std::ios::fixed, std::ios::floatfield );
   }
   vout::fixed_width_output result(res,vout::fw_size);

   // inverse number of data samples * muB
   const double ic = 0.5*( 2.0 /(3.0*constants::kB_eV)) / mean_counter;

   for(int mask_id=0; mask_id<mask_size; ++mask_id){
      if(header){
         result << name + std::to_string(mask_id) + "_mean_Ts";
      }
      else{
         result << mean_lattice_temp[mask_id ]*ic;
      }
   }

   return result.str();

}

} // end of namespace stats
