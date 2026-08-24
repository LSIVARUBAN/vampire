//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Lennon Sivaruban 2026. All rights reserved.
//
//   Email: lennon.sivaruban@postgrad.manchester.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers

// Vampire headers
#include "anisotropy.hpp"

// anisotropy module headers
#include "internal.hpp"

namespace anisotropy {
   namespace internal {

      //---------------------------------------------------------------------------------
      // Function to add fourth order Neel anisotropy valid for collinear spins 
      //       
      // Model derived from Eq 6 of:
      // Nieves, P., Tranchida, J., Arapan, S. and Legut, D., 2021. Spin-lattice model for cubic crystals. Physical Review B, 103(9), p.094437.
      //
      // Hamiltonian:
      //    The Neel energy for a pair of atoms (i,j) is given by:
      //    H_Neel = - (1/2) * sum_neighbors { l(r) * [ (s.e)^2 - 1/3 ]  +  q(r) * [ (s.e)^4 - 6/7*(s.e)^2 + 3/35 ] }
      //
      //    where: x = (s . e)
      //
      // Group terms by powers of x:
      //    H_Neel = - (1/2) * sum_neighbors { [ l(r) - 6/7*q(r) ] * x^2   +   [ q(r) ] * x^4   +   constants }
      //
      //    Define effective coefficients:
      //    L_eff = l(r) - (6/7)*q(r)
      //    Q_eff = q(r)
      //
      // Tensor Definitions (calculated in initialize_neel.cpp):
      //    T2_ab     = sum_neighbors { L_eff * e_a * e_b }             (Rank 2 Tensor a,b,c = x,y,z)
      //    T4_abcd   = sum_neighbors { Q_eff * e_a * e_b * e_c * e_d } (Rank 4 Tensor a,b,c,d = x,y,z)
      //
      //    Substituting tensors into the Hamiltonian:
      //    H_Neel = - (1/2) * [ (s . T2 . s)  +  (s . T4 . s . s . s) ]
      //
      // Field Calculation: H_eff = - dE/ds (energy units so no 1/μ):
      //    A) Quadratic Term:
      //           d/ds(Quadratic Term):
      //               d/ds [ -1/2 * (s . T2 . s) ] 
      //               = -1/2 * d/ds [ sum_ab ( T2_ab * s_a * s_b ) ]
      //               = -1/2 * 2 * (T2 . s) 
      //               = - (T2 . s) 
      //           Therefore H_quad = - dH_quad/ds = (T2 . s)
      //
      //    B) Quartic Term:
      //           d/ds [ -1/2 * (s . T4 . s . s . s) ]
      //           = -1/2 * d/ds [ sum_abcd ( T4_abcd * s_a * s_b * s_c * s_d ) ]
      //           = -1/2 * 4 * (T4 . s . s . s)
      //            = - 2 * (T4 . s . s . s)
      //
      //           Therefore H_quart = - dH_quart/ds = 2 * (T4 . s . s . s)
      //
      //---------------------------------------------------------------------------------

      void neel_fields_fourth(std::vector<double>& spin_array_x,
                            std::vector<double>& spin_array_y,
                            std::vector<double>& spin_array_z,
                            std::vector<int>&    atom_material_array,
                            std::vector<double>& field_array_x,
                            std::vector<double>& field_array_y,
                            std::vector<double>& field_array_z,
                            const int start_index,
                            const int end_index){

         if(!internal::enable_neel_fourth_order_anisotropy) return;

        // loop over all atoms
         for(int atom = start_index; atom < end_index; atom++){

            // Temporary arrays to hold spin vector and field update
            const double s[3] = {spin_array_x[atom],
                                 spin_array_y[atom],
                                 spin_array_z[atom]};

            // Field update is used for every atom so reset it here
            double field_update[3] = {0.0,
                                      0.0,
                                      0.0};

            const int idx2 = 9 * atom;  // Starting index for Rank 2 Tensor
            const int idx4 = 81 * atom; // Starting index for Rank 4 Tensor

            // ---------------------------------------------------------
            // Calculate Quadratic Contribution (Rank 2 Tensor)
            // H_alpha += sum_beta ( T2_ab * s_b ) corresponds to: H = (T2 . s)
            // ---------------------------------------------------------
            for (int a = 0; a < 3; a++) { // alpha (field component index x, y, z)
               double h_component = 0.0;
               
               // h_a = sum_b ( T_ab * s_b )
               for (int b = 0; b < 3; b++) { // beta (spin component index)
                  h_component += internal::neel_tensor[idx2 + 3*a + b] * s[b];
               }
               
               // Add quadratic contribution to local site accumulator
               field_update[a] += h_component; 
            }

            // ---------------------------------------------------------
            // Calculate Quartic Contribution (Rank 4 Tensor)
            // H_alpha += 2.0 * sum_bcd ( T4_abcd * s_b * s_c * s_d ) corresponds to: H = 2 * (T4 . s . s . s)
            // ---------------------------------------------------------
            for (int a = 0; a < 3; a++) { // alpha (field component index x, y, z)
               double h_component = 0.0;
               
               // h_a = sum_b sum_c sum_d ( T_abcd * s_b * s_c * s_d )
               for (int b = 0; b < 3; b++) {
                  for (int c = 0; c < 3; c++) {
                     for (int d = 0; d < 3; d++) {
                        
                        // Map 4D index [a][b][c][d] onto 1D flat array index
                        int flat_k = 27*a + 9*b + 3*c + d;
                        
                        h_component += internal::neel_tensor_4[idx4 + flat_k] * s[b] * s[c] * s[d];
                     }
                  }
               }
               
               // Add quartic contribution to local site accumulator
               field_update[a] += 2.0 * h_component; 
            }

            // Add the total local field to the global field arrays
            field_array_x[atom] += field_update[0];
            field_array_y[atom] += field_update[1];
            field_array_z[atom] += field_update[2];
         }
      }

      //---------------------------------------------------------------------------------
      // Function to add Neel Energy
      // H_Neel = - (1/2) * [ (s . T2 . s)  +  (s . T4 . s . s . s) ]
      //---------------------------------------------------------------------------------
      double neel_energy_fourth(const int atom,
                              const int mat,
                              const double sx,
                              const double sy,
                              const double sz){

         double s[3] = {sx, sy, sz};
         const int idx2 = 9 * atom; // Rank 2 tensor start index
         const int idx4 = 81 * atom; // Rank 4 tensor start index

         double energy_quad = 0.0; // Quadratic energy term
         double energy_quart = 0.0; // Quartic energy term

         // Quadratic Energy Term: sum_ab ( s_a * T2_ab * s_b )
         for(int a=0; a<3; a++){
            for(int b=0; b<3; b++){
               energy_quad += s[a] * internal::neel_tensor[idx2 + 3*a + b] * s[b];
            }
         }

         // Quartic Energy Term: sum_abcd ( s_a * T4_abcd * s_b * s_c * s_d )
         for(int a=0; a<3; a++){
            for(int b=0; b<3; b++){
               for(int c=0; c<3; c++){
                  for(int d=0; d<3; d++){
                     int flat_k = 27*a + 9*b + 3*c + d; // 4D tensor index flattened to 1D
                     energy_quart += s[a] * internal::neel_tensor_4[idx4 + flat_k] * s[b] * s[c] * s[d];
                  }
               }
            }
         }

         // Combine terms, apply the -0.5 prefactor from the Hamiltonian and add the constants
         return -0.5 * (energy_quad + energy_quart) + internal::neel_fourth_order_energy_constant[atom];
      }

   } // end of internal namespace
} // end of anisotropy namespace