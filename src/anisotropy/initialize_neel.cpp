//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sam Westmoreland and Richard Evans 2017. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <map>
#include <cstring>
#include <limits>

// Vampire headers
#include "atoms.hpp" // to be removed
#include "create.hpp" // to be removed
#include "anisotropy.hpp"
#include "errors.hpp"
#include "units.hpp"
#include "vio.hpp"
#include "sim.hpp"  // for sim::constraint_theta, sim::constraint_phi

// anisotropy module headers
#include "internal.hpp"

namespace anisotropy{

namespace internal{

std::string neel_pairs_csv_name; // sim::constraint_theta and sim::constraint_phi as "theta_XXXdeg_phi_YYYdeg.csv"

// Helper to map material id to the labels used in the Python CSVs
static inline const char* mat_label(unsigned int id) {
    switch (id) {
        case 0: return "FeA";
        case 1: return "FeB";
        case 2: return "O";
        default: return "X";
    }
}

// Helper to write per-pair CSV in the current directory (to compare with Python scripts)
// Columns:
// fe_id,o_id,fe_site,fe_coord_z,ex,ey,ez,sx,sy,sz,k_J,S_dot_e,pair_energy_J,fe_x,fe_y,fe_z,o_x,o_y,o_z
static void write_neel_pairs_csv(
    const std::string& csv_name,
    const std::vector<std::vector<bool>>& nnmask,
    const std::vector<std::vector<neighbours::neighbour_t>>& cneighbourlist
) {
    const std::string fname = csv_name.empty() ? std::string("results_pairs.csv") : csv_name;
    std::ofstream out(fname);
    if (!out) {
        zlog << zTs() << "Error opening " << fname << " for writing." << std::endl;
        return;
    }
    out << "fe_id,o_id,fe_site,fe_coord_z,"
           "ex,ey,ez,sx,sy,sz,"
           "k_J,S_dot_e,pair_energy_J,"
           "fe_x,fe_y,fe_z,o_x,o_y,o_z\n";
    out << std::setprecision(16);

    // Iterate surface Fe atoms
    for (int i = 0; i < atoms::num_atoms; ++i) {
        if (!atoms::surface_array[i]) continue;

        const unsigned int imat = atoms::type_array[i];
        const char* ilabel = mat_label(imat);
        if (!(std::strcmp(ilabel, "FeA") == 0 || std::strcmp(ilabel, "FeB") == 0)) continue;

        // Count Fe–O nearest neighbours for coordination z
        unsigned int zcoord = 0;
        for (unsigned int nn = 0; nn < cneighbourlist[i].size(); ++nn) {
            if (!nnmask[i][nn]) continue;
            const unsigned int jz = cneighbourlist[i][nn].nn;
            if (mat_label(atoms::type_array[jz])[0] == 'O') ++zcoord;
        }

        // Spins for atom i
        const double sx = atoms::x_spin_array[i];
        const double sy = atoms::y_spin_array[i];
        const double sz = atoms::z_spin_array[i];

        // Write one row per Fe–O pair
        for (unsigned int nn = 0; nn < cneighbourlist[i].size(); ++nn) {
            if (!nnmask[i][nn]) continue;

            const unsigned int j = cneighbourlist[i][nn].nn;
            const unsigned int jmat = atoms::type_array[j];
            if (mat_label(jmat)[0] != 'O') continue; // Fe–O only

            // Vector i->j and unit e_ij
            double vx = cneighbourlist[i][nn].vx;
            double vy = cneighbourlist[i][nn].vy;
            double vz = cneighbourlist[i][nn].vz;
            const double rij = std::sqrt(vx*vx + vy*vy + vz*vz);
            if (rij <= 0.0) continue;
            const double invrij = 1.0 / rij;
            const double ex = vx * invrij;
            const double ey = vy * invrij;
            const double ez = vz * invrij;

            // k is in Joules assuming 1/2 is taken into account by mp[imat].kij[jmat]
            double k_J = anisotropy::internal::mp[imat].kij[jmat];
            if (internal::neel_range_dependent) {
                // exponential factor if needed
                k_J *= std::exp(-neel_exponential_factor * (rij - neel_exponential_range) / neel_exponential_range);
            }

            const double SdotE = sx*ex + sy*ey + sz*ez;
            const double pair_energy_J = -k_J * (SdotE * SdotE);

            // Positions
            const double ix = atoms::x_coord_array[i];
            const double iy = atoms::y_coord_array[i];
            const double iz = atoms::z_coord_array[i];
            const double jx = atoms::x_coord_array[j];
            const double jy = atoms::y_coord_array[j];
            const double jz = atoms::z_coord_array[j];

            // Row
            out << i << ',' << j << ','
                << ilabel << ',' << zcoord << ','
                << ex << ',' << ey << ',' << ez << ','
                << sx << ',' << sy << ',' << sz << ','
                << k_J << ',' << SdotE << ',' << pair_energy_J << ','
                << ix << ',' << iy << ',' << iz << ','
                << jx << ',' << jy << ',' << jz << '\n';
        }
    }

    out.close();
    zlog << zTs() << "Wrote pairs CSV: " << fname << std::endl;
}

    //---------------------------------------------------------------------------
    // Function to calculate fourth order Neel anisotropy tensors
    //
    // Calculates the per-atom tensors:
    //    T2_ab     = sum_neighbors { L_eff * e_a * e_b }             (Rank 2 Tensor a,b,c = x,y,z)
    //    T4_abcd   = sum_neighbors { Q_eff * e_a * e_b * e_c * e_d } (Rank 4 Tensor a,b,c,d = x,y,z)
    //
    // where:
    //   L_eff = l - (6/7) q   and   Q_eff = q
    //---------------------------------------------------------------------------
    void initialise_neel_fourth_order_anisotropy_tensors(std::vector <std::vector <bool> >& nearest_neighbour_interactions_list,
                                                                          std::vector<std::vector <neighbours::neighbour_t> >& cneighbourlist){

        // Print informative message to log file
        zlog << zTs() << "Initialising fourth order Néel anisotropy tensors for atoms with < threshold number of neighbours." << std::endl;

        // allocate memory for neel anisotropy tensors
        internal::neel_tensor.resize( 9 * atoms::num_atoms, 0.0 );
        internal::neel_tensor_4.resize( 81 * atoms::num_atoms, 0.0 );
        internal::neel_fourth_order_energy_constant.resize( atoms::num_atoms, 0.0 );

        // temporary tensors for calculating sums
        std::vector<double> tmp_tensor2(9, 0.0);
        std::vector<double> tmp_tensor4(81, 0.0);

        //	Populate surface atom and 1D nearest neighbour list and index arrays
        for(int atom = 0; atom < atoms::num_atoms; ++atom){

            // Only calculate parameters for atoms with less than full nn coordination
            if(!atoms::surface_array[atom]) continue;

            // get material id for atom
            const unsigned int imat = atoms::type_array[atom]; 

            // get inverse moment
            const double i_mu_s = 1.0/mp::material[imat].mu_s_SI; 

            // zero tensors
            std::fill(tmp_tensor2.begin(), tmp_tensor2.end(), 0.0);  
            std::fill(tmp_tensor4.begin(), tmp_tensor4.end(), 0.0);

            // accumulate energy offset from constants in Hamiltonian
            double energy_constant = 0.0;

            // loop over all neighbours
            for(unsigned int nn = 0; nn < cneighbourlist[atom].size(); ++nn){

                // only add nearest neighbours to list
                if(nearest_neighbour_interactions_list[atom][nn] != true) continue; 

                // get atom number for neighbour
                const unsigned int natom = cneighbourlist[atom][nn].nn; 

                // get material id for j atom
                const unsigned int jmat = atoms::type_array[natom];

                // Fe–O only in tensor 
                static constexpr unsigned int O_ID = 2; // Hardcoded O material id to be later generalised
                if(jmat != O_ID) continue;

                // get atomic position vector i->j
                double eij[3] = {cneighbourlist[atom][nn].vx, cneighbourlist[atom][nn].vy, cneighbourlist[atom][nn].vz};

                // normalise to unit vector
                const double rij = std::sqrt(eij[0]*eij[0] + eij[1]*eij[1] + eij[2]*eij[2]);
                if(rij <= 0.0) continue;
                const double invrij = 1.0/rij;
                eij[0] = eij[0] * invrij;
                eij[1] = eij[1] * invrij;
                eij[2] = eij[2] * invrij;

                // Convert neel coefficients from J to Tesla
                double l_T = anisotropy::internal::mp[imat].kij[jmat] * i_mu_s;
                double q_T = anisotropy::internal::mp[imat].qij[jmat] * i_mu_s;

                // Assume coefficients have identical exponential distance dependence as simple model if range dependent model is selected
                if(internal::neel_range_dependent){
                    const double factor = std::exp(-neel_exponential_factor * (rij - neel_exponential_range) / neel_exponential_range);
                    l_T *= factor;
                    q_T *= factor;
                }

                // Add the energy_constant = (1/6)*l  - (3/70)*q per neighbour from the Hamiltonian.
                //  H_neel = -(1/2) * sum_j { l[(s.e)^2 - 1/3] + q[(s.e)^4 - 6/7 (s.e)^2 + 3/35] }
                energy_constant += (1.0/6.0) * l_T - (3.0/70.0) * q_T;

                // effective coefficients
                const double L_eff = l_T - (6.0/7.0) * q_T;
                const double Q_eff = q_T;

                // Rank 2 tensor accumulation
                for(int a = 0; a < 3; ++a){
                    for(int b = 0; b < 3; ++b){
                        tmp_tensor2[3*a + b] += L_eff * eij[a] * eij[b];
                    }
                }

                // Rank 4 tensor accumulation
                for(int a = 0; a < 3; ++a){
                    for(int b = 0; b < 3; ++b){
                        for(int c = 0; c < 3; ++c){
                            for(int d = 0; d < 3; ++d){
                                const int flat_k = 27*a + 9*b + 3*c + d;
                                tmp_tensor4[flat_k] += Q_eff * eij[a] * eij[b] * eij[c] * eij[d];
                            }
                        }
                    }
                }

            } // end of neighbour loop

            // Save tensors for atom
            for(int i = 0; i < 9; ++i){
                anisotropy::internal::neel_tensor[9*atom + i] = tmp_tensor2[i];
            }
            for(int i = 0; i < 81; ++i){
                anisotropy::internal::neel_tensor_4[81*atom + i] = tmp_tensor4[i];
            }

            anisotropy::internal::neel_fourth_order_energy_constant[atom] = energy_constant;

        } // end atom loop

        return;
    }

   //---------------------------------------------------------------------------
   // Function to calculate surface anisotropy tensor
   //---------------------------------------------------------------------------
   void initialise_neel_anisotropy_tensor(std::vector <std::vector <bool> >& nearest_neighbour_interactions_list,
                                          std::vector<std::vector <neighbours::neighbour_t> >& cneighbourlist){

      // Print informative message to log file
      zlog << zTs() << "Using Néel pair anisotropy for atoms with < threshold number of neighbours." << std::endl;

      // allocate memory for neel anisotropy tensor
      internal::neel_tensor.resize( 9 * atoms::num_atoms, 0.0 );

      // temporary tensor for calculating sum
      std::vector<double> tmp_tensor(9,0.0);

      //	Populate surface atom and 1D nearest neighbour list and index arrays
      for(int atom=0;atom<atoms::num_atoms;++atom){


         // Only calculate parameters for atoms with less than full nn coordination
         if(atoms::surface_array[atom]){

            // get material id for atom
            const unsigned int imat = atoms::type_array[atom];

            // get inverse moment
            const double i_mu_s = 1.0/mp::material[imat].mu_s_SI;

            // zero tensor
            for(int idx = 0; idx < 9; idx++) tmp_tensor[idx] = 0.0;

            // loop over all neighbours
            for(unsigned int nn = 0; nn < cneighbourlist[atom].size(); nn++){

               // only add nearest neighbours to list
               if(nearest_neighbour_interactions_list[atom][nn]==true){

                  // get atom number for neighbour
                  const unsigned int natom = cneighbourlist[atom][nn].nn;

                  // get material id for j atom
                  const unsigned int jmat = atoms::type_array[natom];

                  // Fe–O only in tensor
                  static constexpr unsigned int O_ID = 2; 
                  if (jmat != O_ID) continue;

                  // get atomic position vector i->j
                  double eij[3]={cneighbourlist[atom][nn].vx, cneighbourlist[atom][nn].vy, cneighbourlist[atom][nn].vz};

                  // normalise to unit vector
                  const double rij = sqrt(eij[0]*eij[0]+eij[1]*eij[1]+eij[2]*eij[2]);
                  const double invrij = 1.0/rij;

                  // normalise components to unit vector
                  eij[0] = eij[0] * invrij;
                  eij[1] = eij[1] * invrij;
                  eij[2] = eij[2] * invrij;

                  // get pair anisotropy constant between atom i and j in Tesla
                  double kij = anisotropy::internal::mp[imat].kij[jmat]*i_mu_s;

                  // Adjust kij by exponential factor if needed
                  if(internal::neel_range_dependent){
                     // Normalised exponential e(r0) = 1.0
                     kij = kij*exp(-neel_exponential_factor*( rij - neel_exponential_range )/neel_exponential_range);
                  };

                  // loop over tensor components and sum total including local anisotropy constant
                  // note inclusion of factor - 2 . 1/2 = -1 in tensor field compared to energy due to derivative
                  for(int i = 0; i < 3; ++i){
                     for(int j = 0; j < 3; ++j){
                        tmp_tensor[ 3*i + j ] += kij * eij[i] * eij[j];
                     }
                  }

                  //int natom = cneighbourlist[atom][nn].nn;
                  //std::cout << "nn_id: " << nn << " j: " << cneighbourlist[atom][nn].nn << "\trange: " << 1.0/invrij << " ";
                  //std::cout << "eij: " << eij[0] << " " << eij[1] << " " << eij[2] << "\tatomi: ";
                  //std::cout << catom_array[atom].x << " " << catom_array[atom].y << " " << catom_array[atom].z << "\tatomj: ";
                  //std::cout << catom_array[natom].x << " " << catom_array[natom].y << " " << catom_array[natom].z << std::endl;
                  //std::cout << atom << "\t" << natom << "\t" << kij << "\t" << anisotropy::internal::mp[imat].kij[jmat] << "\t" << i_mu_s << std::endl;

               }

            } // end of neighbour loop

            // save tensor for atom
            anisotropy::internal::neel_tensor[ 9 * atom + 0 ] = tmp_tensor [0];
            anisotropy::internal::neel_tensor[ 9 * atom + 1 ] = tmp_tensor [1];
            anisotropy::internal::neel_tensor[ 9 * atom + 2 ] = tmp_tensor [2];

            anisotropy::internal::neel_tensor[ 9 * atom + 3 ] = tmp_tensor [3];
            anisotropy::internal::neel_tensor[ 9 * atom + 4 ] = tmp_tensor [4];
            anisotropy::internal::neel_tensor[ 9 * atom + 5 ] = tmp_tensor [5];

            anisotropy::internal::neel_tensor[ 9 * atom + 6 ] = tmp_tensor [6];
            anisotropy::internal::neel_tensor[ 9 * atom + 7 ] = tmp_tensor [7];
            anisotropy::internal::neel_tensor[ 9 * atom + 8 ] = tmp_tensor [8];

         }
         //std::cout << std::endl;
      } // end of atom loop

      //---------------------------------------------------------------------------------
      // Output Neel tensors to file
      //---------------------------------------------------------------------------------

      std::ofstream ofile("Neel_tensor.txt");

      for(int atom=0; atom < atoms::num_atoms; ++atom){
         // only output tensors for surface atoms
         if(atoms::surface_array[atom]){
            ofile << atom << "\t" << atoms::x_coord_array[atom] << "\t" << atoms::y_coord_array[atom] << "\t" << atoms::z_coord_array[atom] << "\t";
            for(int i=0; i<9; i++) ofile << anisotropy::internal::neel_tensor[ 9 * atom + i ] << "\t";
            ofile << std::endl;
         }
      }

      ofile.close();

      // Output the per-pair CSV for comparison with Python scripts
      std::string csv_name_final = neel_pairs_csv_name;

      if (csv_name_final.empty()) {
          double th_deg = 0.0;
          double ph_deg = 0.0;

          th_deg = sim::constraint_theta;
          ph_deg = sim::constraint_phi;

          auto norm360 = [](double a)->double {
              double x = std::fmod(a, 360.0);
              if (x < 0.0) x += 360.0;
              return x;
          };

          const int th_i = static_cast<int>(std::llround(norm360(th_deg)));
          const int ph_i = static_cast<int>(std::llround(norm360(ph_deg)));

          std::ostringstream oss;
          oss << "theta_" << std::setw(3) << std::setfill('0') << th_i
              << "deg_phi_" << std::setw(3) << std::setfill('0') << ph_i
              << "deg.csv";
          csv_name_final = oss.str();
      }

      // Uncomment to write a CSV for this (theta,phi) run in the current directory.
      // write_neel_pairs_csv(csv_name_final, nearest_neighbour_interactions_list, cneighbourlist);

   } // end of surface anisotropy initialisation

} // end of internal namespace

} // end of anisotropy namespace
