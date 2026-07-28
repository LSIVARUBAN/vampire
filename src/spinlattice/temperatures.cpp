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
#include <vector>

// Vampire headers
#include "constants.hpp"
#include "material.hpp"
#include "sld.hpp"
#include "vmpi.hpp"

// sld module headers
#include "internal.hpp"


namespace sld{


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
                double SH=0.0;
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
                     // weight each atom separately so multi material systems use the correct moment rather than the moment of material zero
                     SxH2 += mu_s_array[type_array[at]] * (SxHx*SxHx + SxHy*SxHy + SxHz*SxHz);
                     SH  = SH +  Sx * Hx + Sy * Hy + Sz*Hz;

                }
               #ifdef MPICF
                  double sums[2] = {SxH2, SH};
                  MPI_Allreduce(MPI_IN_PLACE, sums, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                  SxH2 = sums[0];
                  SH = sums[1];
               #endif

               if(SH <= 0.0) return 0.0;
               double T_spin = 0.5 * constants::muB / constants::kB * SxH2 / SH;

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
            double mass_velocity_squared = 0.0;
            double number_of_atoms = 0.0;

            for (int at=start_index;at<end_index;at++){
               const double atom_mass = sld::internal::mp[type_array[at]].mass.get();
               const double vx = velo_array_x[at];
               const double vy = velo_array_y[at];
               const double vz = velo_array_z[at];
               mass += atom_mass;
               momentum_x += atom_mass * vx;
               momentum_y += atom_mass * vy;
               momentum_z += atom_mass * vz;
               mass_velocity_squared += atom_mass * (vx * vx + vy * vy + vz * vz);
               number_of_atoms += 1.0;
            }

            #ifdef MPICF
               double sums[6] = {mass, momentum_x, momentum_y, momentum_z, mass_velocity_squared, number_of_atoms};
               MPI_Allreduce(MPI_IN_PLACE, sums, 6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
               mass = sums[0];
               momentum_x = sums[1];
               momentum_y = sums[2];
               momentum_z = sums[3];
               mass_velocity_squared = sums[4];
               number_of_atoms = sums[5];
            #endif

            if(mass <= 0.0 || number_of_atoms <= 0.0) return 0.0;

            // remove translation of the whole system before converting its kinetic energy to a phonon temperature
            const double centre_of_mass_term = (momentum_x * momentum_x + momentum_y * momentum_y + momentum_z * momentum_z) / mass;
            const double thermal_mass_velocity_squared = std::max(0.0, mass_velocity_squared-centre_of_mass_term);
            const double T_lat = thermal_mass_velocity_squared / (3.0 * number_of_atoms * constants::kB_eV);

   return T_lat;

}//end of lattice_temperature

   } // end of sld namespace
