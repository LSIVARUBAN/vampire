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
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

// Vampire headers
#include "sld.hpp"
#include "create.hpp"
#include "sim.hpp"
#include "errors.hpp"
#include "material.hpp"

// sld module headers
#include "internal.hpp"


namespace sld{

   namespace{

   // if requested output a force table for the first force evaluation containing contributions just from the mechanical potential
   void write_force_debug_table(const int start_index,
                                const int end_index,
                                const std::vector<double>& x_coord_array,
                                const std::vector<double>& y_coord_array,
                                const std::vector<double>& z_coord_array,
                                const std::vector<double>& initial_forces_x,
                                const std::vector<double>& initial_forces_y,
                                const std::vector<double>& initial_forces_z,
                                const std::vector<double>& forces_array_x,
                                const std::vector<double>& forces_array_y,
                                const std::vector<double>& forces_array_z){

      std::ofstream output("force-debug.data");
      if(!output){
         err::zexit("Unable to open force-debug.data");
      }

      output << "# VAMPIRE mechanical-potential force debug from the first force evaluation\n";
      output << "# atom x_A y_A z_A fx_eV_per_A fy_eV_per_A fz_eV_per_A force_magnitude_eV_per_A\n";
      output << std::setprecision(17);

      for(int atom = start_index; atom < end_index; atom++){
         const int local_atom = atom - start_index;
         const double fx = forces_array_x[atom] - initial_forces_x[local_atom];
         const double fy = forces_array_y[atom] - initial_forces_y[local_atom];
         const double fz = forces_array_z[atom] - initial_forces_z[local_atom];
         const double magnitude = std::sqrt(fx*fx + fy*fy + fz*fz);

         output << atom << " " << x_coord_array[atom] << " " << y_coord_array[atom] << " " << z_coord_array[atom] << " " << fx << " " << fy << " " << fz << " " << magnitude << "\n";
      }
   }

   } // end of anonymous namespace

   void compute_forces(const int start_index, // first atom for exchange interactions to be calculated
               const int end_index,
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
               std::vector<double>& potential_eng){
      const bool write_force_debug = sld::output_force_debug && !sld::internal::force_debug_written;
      std::vector<double> initial_forces_x;
      std::vector<double> initial_forces_y;
      std::vector<double> initial_forces_z;
      // save the forces before the potential adds its contribution, we later subtract these from the final forces to isolate the mechanical potential contribution
      if(write_force_debug){
         const int natoms = end_index - start_index;
         initial_forces_x.resize(natoms);
         initial_forces_y.resize(natoms);
         initial_forces_z.resize(natoms);
         for(int atom = start_index; atom < end_index; atom++){
            const int local_atom = atom - start_index;
            initial_forces_x[local_atom] = forces_array_x[atom];
            initial_forces_y[local_atom] = forces_array_y[atom];
            initial_forces_z[local_atom] = forces_array_z[atom];
         }
      }

      switch(sld::internal::lattice_potential){
         case sld::internal::harmonic_lattice_potential:
            // Harmonic and Morse are simple pair potentials
            internal::compute_forces_harmonic(start_index, end_index,
                                              neighbour_list_start_index, neighbour_list_end_index,
                                              type_array, neighbour_list_array,
                                              x0_coord_array, y0_coord_array, z0_coord_array,
                                              x_coord_array, y_coord_array, z_coord_array,
                                              forces_array_x, forces_array_y, forces_array_z, potential_eng);
            break;

         case sld::internal::morse_lattice_potential:
            internal::compute_forces_morse(start_index, end_index,
                                           neighbour_list_start_index, neighbour_list_end_index,
                                           type_array, neighbour_list_array,
                                           x_coord_array, y_coord_array, z_coord_array,
                                           forces_array_x, forces_array_y, forces_array_z, potential_eng);
            break;

         case sld::internal::snap_lattice_potential:
            // SNAP is many body therefore the force on atom i depends on the neighbour density around each central atom (not just a pair)
            internal::compute_forces_snap(start_index, end_index,
                                          neighbour_list_start_index, neighbour_list_end_index,
                                          type_array, neighbour_list_array,
                                          x_coord_array, y_coord_array, z_coord_array,
                                          forces_array_x, forces_array_y, forces_array_z, potential_eng);
            break;

         case sld::internal::snap_zbl_lattice_potential:
            // SNAP plus a ZBL potential overlay for short range pair repulsion
            internal::compute_forces_snap_zbl(start_index, end_index,
                                             neighbour_list_start_index, neighbour_list_end_index,
                                             type_array, neighbour_list_array,
                                             x_coord_array, y_coord_array, z_coord_array,
                                             forces_array_x, forces_array_y, forces_array_z, potential_eng);
            break;

         default:
            break;
      }

      // write forces from mechnical potential to file once, if requested
      if(write_force_debug){
         write_force_debug_table(start_index, end_index,
                                 x_coord_array, y_coord_array, z_coord_array,
                                 initial_forces_x, initial_forces_y, initial_forces_z,
                                 forces_array_x, forces_array_y, forces_array_z);
         sld::internal::force_debug_written = true;
      }
      
      // calculate THz pulse force
      if(sld::internal::linear_pump_enabled){

        internal::compute_thz(start_index, end_index,
                                x_coord_array, y_coord_array, z_coord_array,
                                forces_array_x, forces_array_y, forces_array_z);
      
        }

      return;

}

namespace internal{

namespace{

struct zbl_coeff_t{
   double d1a; // first exponential decay constant divided by screening length
   double d2a; // second
   double d3a; // third
   double d4a; // fourth
   double zze; // Zi*Zj*e^2/(4*pi*epsilon0) in eV Angstrom
   double sw1; // force switch coefficient A
   double sw2; // force switch coefficient B
   double sw3; // energy switch coefficient A/3
   double sw4; // energy switch coefficient B/4
   double sw5; // energy constant offset
};

double e_zbl(const double r, const zbl_coeff_t& coeff){

   // unswitched screened Coulomb pair energy in eV
   const double c1 = 0.02817;
   const double c2 = 0.28022;
   const double c3 = 0.50986;
   const double c4 = 0.18175;
   const double d1aij = coeff.d1a;
   const double d2aij = coeff.d2a;
   const double d3aij = coeff.d3a;
   const double d4aij = coeff.d4a;
   const double zzeij = coeff.zze;
   const double rinv = 1.0 / r;

   // phi(r/a) is a weighted sum of four exponentials with each term describing one part of the ZBL screening cloud
   double sum = c1 * std::exp(-d1aij * r);
   sum += c2 * std::exp(-d2aij * r);
   sum += c3 * std::exp(-d3aij * r);
   sum += c4 * std::exp(-d4aij * r);

   // E_ZBL(r) = Zi*Zj*e^2/(4*pi*epsilon0) * phi(r/a) / r (in eV)
   double result = zzeij * sum * rinv;

   return result;

}

double dzbldr(const double r, const zbl_coeff_t& coeff){

   // dE_ZBL/dr before the switching derivative is added
   const double c1 = 0.02817;
   const double c2 = 0.28022;
   const double c3 = 0.50986;
   const double c4 = 0.18175;
   const double d1aij = coeff.d1a;
   const double d2aij = coeff.d2a;
   const double d3aij = coeff.d3a;
   const double d4aij = coeff.d4a;
   const double zzeij = coeff.zze;
   const double rinv = 1.0 / r;

   const double e1 = std::exp(-d1aij * r);
   const double e2 = std::exp(-d2aij * r);
   const double e3 = std::exp(-d3aij * r);
   const double e4 = std::exp(-d4aij * r);

   // phi(r/a)
   double sum = c1 * e1;
   sum += c2 * e2;
   sum += c3 * e3;
   sum += c4 * e4;

   // d(phi)/dr
   double sum_p = -c1 * d1aij * e1;
   sum_p -= c2 * d2aij * e2;
   sum_p -= c3 * d3aij * e3;
   sum_p -= c4 * d4aij * e4;

   // dE/dr = zze * (phi'/r - phi/r^2).
   double result = zzeij * (sum_p - sum * rinv) * rinv;

   return result;

}

double d2zbldr2(const double r, const zbl_coeff_t& coeff){

   const double c1 = 0.02817;
   const double c2 = 0.28022;
   const double c3 = 0.50986;
   const double c4 = 0.18175;
   const double d1aij = coeff.d1a;
   const double d2aij = coeff.d2a;
   const double d3aij = coeff.d3a;
   const double d4aij = coeff.d4a;
   const double zzeij = coeff.zze;
   const double rinv = 1.0 / r;

   const double e1 = std::exp(-d1aij * r);
   const double e2 = std::exp(-d2aij * r);
   const double e3 = std::exp(-d3aij * r);
   const double e4 = std::exp(-d4aij * r);

   // phi(r/a)
   double sum = c1 * e1;
   sum += c2 * e2;
   sum += c3 * e3;
   sum += c4 * e4;

   // -d(phi)/dr
   double sum_p = c1 * e1 * d1aij;
   sum_p += c2 * e2 * d2aij;
   sum_p += c3 * e3 * d3aij;
   sum_p += c4 * e4 * d4aij;

   // d2(phi)/dr2
   double sum_pp = c1 * e1 * d1aij * d1aij;
   sum_pp += c2 * e2 * d2aij * d2aij;
   sum_pp += c3 * e3 * d3aij * d3aij;
   sum_pp += c4 * e4 * d4aij * d4aij;

   // d2E/dr2 for E = zze * phi/r.
   double result = zzeij * (sum_pp + 2.0 * sum_p * rinv + 2.0 * sum * rinv * rinv) * rinv;

   return result;

}

zbl_coeff_t set_zbl_coeff(const double zi,
                          const double zj,
                          const double cut_inner,
                          const double cut_outer){

   const double pzbl = 0.23;
   const double a0 = 0.46850;
   const double d1 = 0.20162;
   const double d2 = 0.40290;
   const double d3 = 0.94229;
   const double d4 = 3.19980;
   const double qqr2e_metal = 14.399645;

   zbl_coeff_t coeff;
   // ZBL screening length: a = a0 / (Zi^0.23 + Zj^0.23)
   const double ainv = (std::pow(zi, pzbl) + std::pow(zj, pzbl)) / a0;
   coeff.d1a = d1 * ainv;
   coeff.d2a = d2 * ainv;
   coeff.d3a = d3 * ainv;
   coeff.d4a = d4 * ainv;
   // Coulomb prefactor
   coeff.zze = zi * zj * qqr2e_metal;

   // e      = t^3 * (sw3 + sw4*t) + sw5
   // dedr   = t^2 * (sw1 + sw2*t)
   // d2edr2 = 2*sw1*t + 3*sw2*t^2
   const double tc = cut_outer - cut_inner;
   const double fc = e_zbl(cut_outer, coeff);
   const double fcp = dzbldr(cut_outer, coeff);
   const double fcpp = d2zbldr2(cut_outer, coeff);

   const double swa = (-3.0 * fcp + tc * fcpp) / (tc * tc);
   const double swb = (2.0 * fcp - tc * fcpp) / (tc * tc * tc);
   const double swc = -fc + (tc / 2.0) * fcp - (tc * tc / 12.0) * fcpp;

   coeff.sw1 = swa;
   coeff.sw2 = swb;
   coeff.sw3 = swa / 3.0;
   coeff.sw4 = swb / 4.0;
   coeff.sw5 = swc;

   return coeff;

}

double zbl_single_energy(const double r,
                         const double cut_inner,
                         const zbl_coeff_t& coeff){

   double phi = e_zbl(r, coeff);
   phi += coeff.sw5;

   if(r > cut_inner){
      const double t = r - cut_inner;
      const double eswitch = t*t*t * (coeff.sw3 + coeff.sw4 * t);
      phi += eswitch;
   }

   return phi;

}

// return dE/dr (not force)
double zbl_single_dedr(const double r,
                       const double cut_inner,
                       const zbl_coeff_t& coeff){

   double fforce = dzbldr(r, coeff);

   if(r > cut_inner){
      const double t = r - cut_inner;
      const double fswitch = t*t * (coeff.sw1 + coeff.sw2 * t);
      fforce += fswitch;
   }

   return fforce;

}

std::vector<int> count_neighbours_within_cutoff(const int start_index,
                                                const int end_index,
                                                const std::vector<int>& neighbour_list_start_index,
                                                const std::vector<int>& neighbour_list_end_index,
                                                const std::vector<int>& neighbour_list_array,
                                                const std::vector<double>& x_coord_array,
                                                const std::vector<double>& y_coord_array,
                                                const std::vector<double>& z_coord_array,
                                                const double cutoff){

   std::vector<int> neighbour_counts(end_index - start_index, 0);
   const double cutoff_sqr = cutoff * cutoff;

   for(int i = start_index; i < end_index; i++){
      const double rx = x_coord_array[i];
      const double ry = y_coord_array[i];
      const double rz = z_coord_array[i];
      const int nbr_start = neighbour_list_start_index[i];
      const int nbr_end = neighbour_list_end_index[i] + 1;

      for(int n = nbr_start; n < nbr_end; n++){
         const int j = neighbour_list_array[n];
         if(j == i) continue;

         double dx = rx - x_coord_array[j];
         double dy = ry - y_coord_array[j];
         double dz = rz - z_coord_array[j];

         dx = sld::PBC_wrap(dx, cs::system_dimensions[0], cs::pbc[0]);
         dy = sld::PBC_wrap(dy, cs::system_dimensions[1], cs::pbc[1]);
         dz = sld::PBC_wrap(dz, cs::system_dimensions[2], cs::pbc[2]);

         const double rsq = dx*dx + dy*dy + dz*dz;
         if(rsq > 0.0 && rsq < cutoff_sqr) neighbour_counts[i - start_index]++;
      }
   }

   return neighbour_counts;

}

} // end of anonymous namespace

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
                               const std::vector<double>& forces_array_z){

   const int natoms = end_index - start_index;
   if(natoms <= 0) return;

   int min_neigh = std::numeric_limits<int>::max();
   int max_neigh = 0;
   double avg_neigh = 0.0;
   double min_force = std::numeric_limits<double>::max();
   double max_force = 0.0;
   double avg_force = 0.0;
   int nonfinite_force = 0;

   for(int atom = start_index; atom < end_index; atom++){
      const int local_atom = atom - start_index;
      const int ninside = neighbour_counts[local_atom];
      min_neigh = std::min(min_neigh, ninside);
      max_neigh = std::max(max_neigh, ninside);
      avg_neigh += double(ninside);

      const double fx = forces_array_x[atom] - initial_forces_x[local_atom];
      const double fy = forces_array_y[atom] - initial_forces_y[local_atom];
      const double fz = forces_array_z[atom] - initial_forces_z[local_atom];
      const double force = std::sqrt(fx*fx + fy*fy + fz*fz);

      if(std::isfinite(force)){
         min_force = std::min(min_force, force);
         max_force = std::max(max_force, force);
         avg_force += force;
      }
      else nonfinite_force++;
   }

   avg_neigh /= double(natoms);
   const int finite_force_count = natoms - nonfinite_force;
   if(finite_force_count > 0) avg_force /= double(finite_force_count);
   else min_force = 0.0;

   const std::streamsize old_precision = std::cout.precision();
   std::cout<<std::setprecision(12);
   std::cout<<label<<" force debug call "<<call<<std::endl;
   std::cout<<"neighbours avg/min/max [count]: "<<avg_neigh<<" "<<min_neigh<<" "<<max_neigh<<std::endl;
   std::cout<<"force avg/min/max [eV/Angstrom]: "<<avg_force<<" "<<min_force<<" "<<max_force<<std::endl;
   if(nonfinite_force > 0){
      std::cout<<"nonfinite force count: "<<nonfinite_force<<std::endl;
   }

   std::vector<int> sample_atoms;
   sample_atoms.push_back(start_index);
   if(natoms > 2) sample_atoms.push_back(start_index + natoms / 2);
   if(natoms > 1) sample_atoms.push_back(end_index - 1);

   for(size_t isample = 0; isample < sample_atoms.size(); isample++){
      const int atom = sample_atoms[isample];
      const int local_atom = atom - start_index;
      const double fx = forces_array_x[atom] - initial_forces_x[local_atom];
      const double fy = forces_array_y[atom] - initial_forces_y[local_atom];
      const double fz = forces_array_z[atom] - initial_forces_z[local_atom];
      const double force = std::sqrt(fx*fx + fy*fy + fz*fz);

      std::cout<<"atom "<<atom
               <<" neighbours [count]: "<<neighbour_counts[local_atom]
               <<" force [eV/Angstrom]: "
               <<fx<<" "<<fy<<" "<<fz
               <<" magnitude "<<force<<std::endl;
   }

   std::cout.precision(old_precision);

}


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
            std::vector<double>& potential_eng){



            double rx, ry, rz;
            double rx0, ry0, rz0;
            double dx, dy, dz;
            double dx0, dy0, dz0;
            double fx = 0.0, fy = 0.0, fz = 0.0;
            double rji_sqr, rji, rji0, inv_rji;
            int j, total_int;
            double r_sqr_cut=sld::internal::r_cut_pot*sld::internal::r_cut_pot;
            double energy;
            const int natoms = end_index - start_index;
            const bool debug_this_call = sld::internal::harmonic_debug_enabled &&
                                         sld::internal::harmonic_debug_force_calls <
                                         sld::internal::harmonic_debug_max_force_calls;
            std::vector<int> neighbour_counts;
            std::vector<double> initial_forces_x;
            std::vector<double> initial_forces_y;
            std::vector<double> initial_forces_z;

            if(debug_this_call){
               neighbour_counts.assign(natoms, 0);
               initial_forces_x.assign(natoms, 0.0);
               initial_forces_y.assign(natoms, 0.0);
               initial_forces_z.assign(natoms, 0.0);

               for(int i = start_index; i < end_index; i++){
                  const int local_atom = i - start_index;
                  initial_forces_x[local_atom] = forces_array_x[i];
                  initial_forces_y[local_atom] = forces_array_y[i];
                  initial_forces_z[local_atom] = forces_array_z[i];
               }
            }



            for(int i=start_index;i< end_index; ++i){

               fx = 0.0;
               fy = 0.0;
               fz = 0.0;
               energy=0.0;
               total_int=0;

               rx = x_coord_array[i];
               ry = y_coord_array[i];
               rz = z_coord_array[i];
               rx0 = x0_coord_array[i];
               ry0 = y0_coord_array[i];
               rz0 = z0_coord_array[i];

                //note for sld_neighbour_list_array
                // int nbr_end = neighbour_list_end_index[i];
                //for(int i=start_index;i<= end_index; ++i)
                //for( int n = nbr_start; n <=nbr_end; ++n)



        	      int nbr_start = neighbour_list_start_index[i];
        	      int nbr_end = neighbour_list_end_index[i]+1;

        	      for( int n = nbr_start; n < nbr_end; ++n){
        	        j = neighbour_list_array[n];

        	        if ( j != i){
        		       dx = -x_coord_array[j] + rx;
        		       dy = -y_coord_array[j] + ry;
        		       dz = -z_coord_array[j] + rz;
                   dx0 = -x0_coord_array[j] + rx0;
                   dy0 = -y0_coord_array[j] + ry0;
                   dz0 = -z0_coord_array[j] + rz0;

                   dx = sld::PBC_wrap( dx, cs::system_dimensions[0], cs::pbc[0]);
                   dy = sld::PBC_wrap( dy, cs::system_dimensions[1], cs::pbc[1]);
                   dz = sld::PBC_wrap( dz, cs::system_dimensions[2], cs::pbc[2]);
                   dx0 = sld::PBC_wrap( dx0, cs::system_dimensions[0], cs::pbc[0]);
                   dy0 = sld::PBC_wrap( dy0, cs::system_dimensions[1], cs::pbc[1]);
                   dz0 = sld::PBC_wrap( dz0, cs::system_dimensions[2], cs::pbc[2]);



        		       rji_sqr = dx*dx + dy*dy + dz*dz;

        		       if( rji_sqr < r_sqr_cut){

                       total_int++;


        		           rji = sqrt(rji_sqr);
                           rji0 = sqrt(dx0*dx0 + dy0*dy0 + dz0*dz0);
        		           inv_rji = 1.0/ rji;

                       energy += (rji-rji0)*(rji-rji0);


		               fx -=  (rji-rji0)*dx*inv_rji ; //2 (rji-rj0)*dx*inv_rji -> 2 went at the end
                       fy -=  (rji-rji0)*dy*inv_rji ;
                       fz -=  (rji-rji0)*dz*inv_rji ;


                     }
        	   	}
        	    }
             const unsigned int imat = type_array[i];

             double V0=sld::internal::mp[imat].V0.get(); //0.15

        	    forces_array_x[i] += V0 * 2.0 * fx;
        	    forces_array_y[i] += V0  * 2.0 * fy;
        	    forces_array_z[i] += V0  * 2.0 * fz;
                potential_eng[i] = 0.5 * V0 * energy;

                if(debug_this_call){
                   const int local_atom = i - start_index;
                   neighbour_counts[local_atom] = total_int;
                }

  }

            if(debug_this_call){
               sld::internal::print_force_debug_summary("Harmonic",
                                                        sld::internal::harmonic_debug_force_calls + 1,
                                                        start_index, end_index,
                                                        neighbour_counts,
                                                        initial_forces_x,
                                                        initial_forces_y,
                                                        initial_forces_z,
                                                        forces_array_x,
                                                        forces_array_y,
                                                        forces_array_z);
               sld::internal::harmonic_debug_force_calls++;
            }


     return;
            }

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
            std::vector<double>& potential_eng){



            double rx, ry, rz;
            double dx, dy, dz;
            double fx = 0.0, fy = 0.0, fz = 0.0;
            double rji_sqr, rji, inv_rji; // rji0,
            int j; //, total_int;
            double r_sqr_cut=sld::internal::r_cut_pot*sld::internal::r_cut_pot;
            double energy;

           double alpha_m= sld::internal::alpha_m; //1.3885;
           //double r0_m= sld::internal::r0_m;//2.845;
           double morse_D=sld::internal::morse_D; //0.4174;;
           double morse_beta=sld::internal::morse_beta;//exp( alpha_m * r0_m);
           double morse_factor = sld::internal::morse_factor; //-2.0 * morse_D * alpha_m;




            for(int i=start_index;i< end_index; ++i){

               fx = 0.0;
               fy = 0.0;
               fz = 0.0;
               energy=0.0;
               //total_int=0;


               rx = x_coord_array[i];
               ry = y_coord_array[i];
               rz = z_coord_array[i];


        	      int nbr_start = neighbour_list_start_index[i];
        	      int nbr_end = neighbour_list_end_index[i]+1;

        	      for( int n = nbr_start; n < nbr_end; ++n){
        	        j = neighbour_list_array[n];

        	        if ( j != i){
        		       dx = x_coord_array[j] -rx;
        		       dy = y_coord_array[j]- ry;
        		       dz = z_coord_array[j]- rz;

                   dx = sld::PBC_wrap( dx, cs::system_dimensions[0], cs::pbc[0]);
                   dy = sld::PBC_wrap( dy, cs::system_dimensions[1], cs::pbc[1]);
                   dz = sld::PBC_wrap( dz, cs::system_dimensions[2], cs::pbc[2]);



        		       rji_sqr = dx*dx + dy*dy + dz*dz;

        		       if( rji_sqr < r_sqr_cut){

        		           rji = sqrt(rji_sqr);
    		               inv_rji = 1.0/ rji;

		                   double y = morse_beta * exp( - alpha_m * rji);
  		                   double f_morse = y * ( y - 1.0);

        		           //std::cout<<"embedded "<<i<<"\t"<<j<<"\t"<<rji<<"\t"<<position_rij<<"\t"<<int(position_rij)<<"\t"<<v_rij<<"\t"<<rho_rij<<std::endl;



		                 fx +=  f_morse*dx*inv_rji;
                         fy +=  f_morse*dy*inv_rji;
                         fz +=  f_morse*dz*inv_rji;

                       //if(i==0) std::cout<<"fxyz "<<i <<"\t"<<j<<"\t"<<rji<<"\t"<<f_morse*dx*inv_rji<<"\t"<<f_morse*dy*inv_rji<<"\t"<<f_morse*dz*inv_rji<< std::endl;//<<(rji-rji0)*dx*inv_rji<<"\t"<<(rji-rji0)*dy*inv_rji<<"\t"<<(rji-rji0)*dz*inv_rji<<std::endl;
                       //if(i==0) std::cout<<"fxyz "<<i <<"\t"<<j<<"\t"<<rji<<"\t"<<x_coord_array[j]<<"\t"<<y_coord_array[j]<<"\t"<<z_coord_array[j]<< "\t"<<dx<<"\t"<<dy<<"\t"<<dz<<"\t"<<f_morse<<std::endl;//<<(rji-rji0)*dx*inv_rji<<"\t"<<(rji-rji0)*dy*inv_rji<<"\t"<<(rji-rji0)*dz*inv_rji<<std::endl;

                       //if(i==100) std::cout<<"fxyz "<<i <<"\t"<<j<<"\t"<<type_array[i]<<"\t"<<inv_rji<<std::endl;//<<(rji-rji0)*dx*inv_rji<<"\t"<<(rji-rji0)*dy*inv_rji<<"\t"<<(rji-rji0)*dz*inv_rji<<std::endl;
                         energy += y * ( y - 2.0);

                     }


        	   	}


        	    }





        	    forces_array_x[i] += fx*morse_factor;
        	    forces_array_y[i] += fy*morse_factor;
        	    forces_array_z[i] += fz*morse_factor;

              potential_eng[i] = morse_D * energy;

  }



   return;
          }


// Author: Muhammad Hamza Asim
// Applies a time-dependent terahertz (THz) phonon-driven force to atoms
void compute_thz(const int start_index,
                 const int end_index,
                 const std::vector<double>& x_coord_array,
                 const std::vector<double>& y_coord_array,
                 const std::vector<double>& z_coord_array,
                 std::vector<double>& forces_array_x,
                 std::vector<double>& forces_array_y,
                 std::vector<double>& forces_array_z)
{


    // Calculate the Current Physical Time
    // Converts the current step number into the actual time in seconds.
    uint64_t current_step = sim::time; 
    double dt = mp::dt_SI;             // Get the duration of one step in seconds 
    double current_time = static_cast<double>(current_step) * dt; // Total elapsed time (seconds)


    if (current_time < sld::internal::phonon_pulse_start_time || current_time > sld::internal::phonon_pulse_end_time) {
        return; // The pulse is off, so do nothing.
    }
 
    // Calculate the Raw Force for Each Atom
    // Loop through each atom to calculate the
    // force from the THz pulse and add it to a running total for later averaging.
    double sumx = 0.0, sumy = 0.0, sumz = 0.0;
    std::vector<double> f_thz_temp_x(end_index - start_index); // Temporary storage
    std::vector<double> f_thz_temp_y(end_index - start_index);
    std::vector<double> f_thz_temp_z(end_index - start_index);

    // Pre-calculate angular frequency (ω = 2πν) for efficiency
    const double twopi_niu = sld::internal::phonon_frequency * 6.28318530718;
    for (int i = start_index; i < end_index; ++i) {


        // Calculate the position-dependent part of the wave (k·r)
        double kr = sld::internal::phonon_wavevector[0] * x_coord_array[i] +
                    sld::internal::phonon_wavevector[1] * y_coord_array[i] +
                    sld::internal::phonon_wavevector[2] * z_coord_array[i];

        // Calculate the argument for the main physics equation: cos(ωt - k·r)
        double arg = current_time * twopi_niu - kr;
        double cos_factor = std::cos(arg); // The oscillating part of the force
        double sin_factor = std::sin(arg);
        
        // Calculate the raw force and store it temporarily
        int local_index = i - start_index;
        f_thz_temp_x[local_index] = sld::internal::phonon_force_amplitude[0] * cos_factor;
        f_thz_temp_y[local_index] = sld::internal::phonon_force_amplitude[1] * sin_factor;
        f_thz_temp_z[local_index] = sld::internal::phonon_force_amplitude[2];

        // Add to the running total
        sumx += f_thz_temp_x[local_index];
        sumy += f_thz_temp_y[local_index];
        sumz += f_thz_temp_z[local_index];

        
      }

    // Calculate the Center-of-Mass Correction
    // Calculate the average force and will subtract it from each atom's individual force.
    int counter = end_index - start_index;
    double avg_fx = 0.0, avg_fy = 0.0, avg_fz = 0.0;
    if (counter > 0) {
        avg_fx = sumx / static_cast<double>(counter);
        avg_fy = sumy / static_cast<double>(counter);
        avg_fz = sumz / static_cast<double>(counter);
    }
  
    
    // Apply the Final, Corrected Force 
    // Loop through the atoms again, retrieve the
    // temporarily stored raw force, subtract the average
    for (int i = start_index; i < end_index; ++i) {
        int local_index = i - start_index;
        double corrected_fx = f_thz_temp_x[local_index] - avg_fx;
        double corrected_fy = f_thz_temp_y[local_index] - avg_fy;
        double corrected_fz = f_thz_temp_z[local_index] - avg_fz;

        // Add the final calculated force to the atom's total force
        forces_array_x[i] += corrected_fx;
        forces_array_y[i] += corrected_fy;
        forces_array_z[i] += corrected_fz;

        if (i == 0 && current_step % 200 == 0) {
        std::cout << "  Corrected Force X : " << corrected_fx << " N" << std::endl;
       
      }
    }

    return;
}

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
            std::vector<double>& potential_eng){

   sld::internal::snap_potential.compute_forces(start_index, end_index,
                                                neighbour_list_start_index, neighbour_list_end_index,
                                                type_array, neighbour_list_array,
                                                x_coord_array, y_coord_array, z_coord_array,
                                                forces_array_x, forces_array_y, forces_array_z,
                                                potential_eng);

   return;
}

void compute_forces_zbl_overlay(const int start_index,
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
            std::vector<double>& potential_eng){

   // cut_inner is where the smooth switch begins, cut_outer is where the ZBL contribution has gone to zero
   const double cut_inner = sld::internal::zbl_inner_cutoff;
   const double cut_outer = sld::internal::zbl_outer_cutoff;
   const double cut_outer_sqr = cut_outer * cut_outer;
   const double z = sld::internal::zbl_atomic_number;
   const zbl_coeff_t coeff = set_zbl_coeff(z, z, cut_inner, cut_outer);

   // loop over central atoms owned by this force call
   // the arrays contain field or SNAP contributions so accumulate with the zbl forces
   for(int i = start_index; i < end_index; i++){
      const double rx = x_coord_array[i];
      const double ry = y_coord_array[i];
      const double rz = z_coord_array[i];

      // accumulate only the ZBL force contribution for this central atom
      double fx = 0.0;
      double fy = 0.0;
      double fz = 0.0;
      double energy = 0.0;

      // use the spin-lattice neighbour list
      const int nbr_start = neighbour_list_start_index[i];
      const int nbr_end = neighbour_list_end_index[i] + 1;

      for(int n = nbr_start; n < nbr_end; n++){
         const int j = neighbour_list_array[n];

         if(j != i){
            // vector from neighbour j to central atom i
            double dx = rx - x_coord_array[j];
            double dy = ry - y_coord_array[j];
            double dz = rz - z_coord_array[j];

            dx = sld::PBC_wrap(dx, cs::system_dimensions[0], cs::pbc[0]);
            dy = sld::PBC_wrap(dy, cs::system_dimensions[1], cs::pbc[1]);
            dz = sld::PBC_wrap(dz, cs::system_dimensions[2], cs::pbc[2]);

            const double rsq = dx*dx + dy*dy + dz*dz;

            // only pairs inside the ZBL outer cutoff contribute
            if(rsq > 0.0 && rsq < cut_outer_sqr){
               const double r = std::sqrt(rsq);

               // fpair = -(dE/dr)/r, then multiply by delx, dely, delz to get the force
               const double fpair = -zbl_single_dedr(r, cut_inner, coeff) / r;

               // Cartesian force on central atom i from neighbour j
               fx += dx * fpair;
               fy += dy * fpair;
               fz += dz * fpair;

               // each pair is encountered from both atoms so half the force on each central atom
               energy += 0.5 * zbl_single_energy(r, cut_inner, coeff);
            }
         }
      }

      // overlay the ZBL contribution onto any previous force terms
      forces_array_x[i] += fx;
      forces_array_y[i] += fy;
      forces_array_z[i] += fz;
      potential_eng[i] += energy;
   }

}

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
            std::vector<double>& potential_eng){

   const int natoms = end_index - start_index;
   const bool debug_this_call = sld::internal::snap_potential.debug_this_call();
   std::vector<int> neighbour_counts;
   std::vector<double> initial_forces_x;
   std::vector<double> initial_forces_y;
   std::vector<double> initial_forces_z;

   if(debug_this_call){
      neighbour_counts = count_neighbours_within_cutoff(start_index, end_index,
                                                        neighbour_list_start_index,
                                                        neighbour_list_end_index,
                                                        neighbour_list_array,
                                                        x_coord_array,
                                                        y_coord_array,
                                                        z_coord_array,
                                                        sld::internal::r_cut_pot);
      initial_forces_x.assign(natoms, 0.0);
      initial_forces_y.assign(natoms, 0.0);
      initial_forces_z.assign(natoms, 0.0);

      for(int i = start_index; i < end_index; i++){
         const int local_atom = i - start_index;
         initial_forces_x[local_atom] = forces_array_x[i];
         initial_forces_y[local_atom] = forces_array_y[i];
         initial_forces_z[local_atom] = forces_array_z[i];
      }
   }

   // calculate SNAP and fill potential_eng with SNAP atomic energies and add SNAP forces
   sld::internal::snap_potential.compute_forces(start_index, end_index,
                                                neighbour_list_start_index, neighbour_list_end_index,
                                                type_array, neighbour_list_array,
                                                x_coord_array, y_coord_array, z_coord_array,
                                                forces_array_x, forces_array_y, forces_array_z,
                                                potential_eng,
                                                false);

   // add the pairwise ZBL repulsion
   sld::internal::compute_forces_zbl_overlay(start_index, end_index,
                                             neighbour_list_start_index, neighbour_list_end_index,
                                             type_array, neighbour_list_array,
                                             x_coord_array, y_coord_array, z_coord_array,
                                             forces_array_x, forces_array_y, forces_array_z,
                                             potential_eng);

   if(debug_this_call){
      sld::internal::print_force_debug_summary("SNAP+ZBL",
                                               sld::internal::snap_potential.debug_call_number(),
                                               start_index, end_index,
                                               neighbour_counts,
                                               initial_forces_x,
                                               initial_forces_y,
                                               initial_forces_z,
                                               forces_array_x,
                                               forces_array_y,
                                               forces_array_z);
      sld::internal::snap_potential.increment_debug_force_calls();
   }

}

         } //end of internal
      } // end of sld namespace
