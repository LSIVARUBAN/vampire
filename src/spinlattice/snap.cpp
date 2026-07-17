//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   Dev: Lennon Sivaruban 2026
//
//   Email: lennon.sivaruban@postgrad.manchester.ac.uk
//
//   SNAP descriptor implementation adapted from the LAMMPS SNA implementation
//
//------------------------------------------------------------------------------

/*
To be added:
Different neighbour list for SNAP, rather than using the Vampire neighbour list so that the SNAP cutoff can be larger than the Vampire cutoff.
Multi element support
Parallel integrator support for SNAP/MLIAPs
Quadratic SNAP support (improves accuracy?)
*/

// standard library headers
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>

// vampire headers
#include "create.hpp"
#include "errors.hpp"
#include "sld.hpp"
#include "vio.hpp"

// sld headers
#include "internal.hpp"
#include "snap.hpp"

namespace{

   const double pi = 3.14159265359;

   // SNAP potential files follow the LAMMPS file format
   // These helpers are used to parse the snap input files

   // trim whitespace from the start and end of a string
   std::string trim(const std::string& text){

      const std::string whitespace = " \t\n\r\f";
      const std::string::size_type first = text.find_first_not_of(whitespace);
      if(first == std::string::npos) return "";
      const std::string::size_type last = text.find_last_not_of(whitespace);
      return text.substr(first, last - first + 1);

   }

   // trim comments from a string
   std::string trim_comment(const std::string& text){

      return trim(text.substr(0, text.find('#')));

   }

   // get words from a line, ignoring comments and whitespace
   std::vector<std::string> words_from_line(const std::string& line){

      std::vector<std::string> words;
      std::istringstream source(trim_comment(line));
      std::string word;
      while(source >> word) words.push_back(word);
      return words;

   }

   // read words from a file, ignoring comments and whitespace
   std::vector<std::string> read_words(std::ifstream& file){

      std::string line;
      while(std::getline(file, line)){
         std::vector<std::string> words = words_from_line(line);
         if(words.size() != 0) return words;
      }

      return std::vector<std::string>();

   }

   // make sure a string is formatted as a numeric value
   bool valid_numeric_token(const std::string& word){

      std::istringstream source(word);
      double value = 0.0;
      source >> value;
      return source && source.eof();

   }

   // numeric conversion from string to double
   double numeric(const std::string& word, const std::string& filename){

      if(!valid_numeric_token(word)){
         err::zexit("Invalid numeric value \"" + word + "\" in SNAP file \"" + filename + "\"");
      }
      return vin::str_to_double(word);

   }


   // numeric conversion from string to int
   int inumeric(const std::string& word, const std::string& filename){

      std::istringstream source(word);
      int value = 0;
      source >> value;
      if(!source || !source.eof()){
         err::zexit("Invalid integer value \"" + word + "\" in SNAP file \"" + filename + "\"");
      }
      return value;

   }

} // end of anonymous namespace

namespace sld{

   namespace internal{

      snap_sna_t::snap_sna_t() : 
         ncoeff(0), // number of bispectrum coefficients (excluding the coefficient for j1=j2=j=0)
         nmax(0), 
         twojmax(0), // larger values create more B components
         rfac0(0.99363), // angular scaling factor for 3D-4D mapping
         rmin0(0.0), // inner radius for the SNAP cutoff and coordinate mapping
         switch_flag(1), // smooth neighbour cutoff
         bzero_flag(1), // subtracts isolated-atom bispectrum offsets
         idxcg_max(0), // maximum index for the Clebsch-Gordan coefficient list
         idxu_max(0), // maximum index for the U expansion list
         idxz_max(0), // maximum index for the Z intermediate list
         idxb_max(0) // maximum index for the B bispectrum list
      {
      }

      // initialise the bispectrum descriptor
      void snap_sna_t::initialise(const double rfac0_in,
                                  const int twojmax_in,
                                  const double rmin0_in,
                                  const int switch_flag_in,
                                  const int bzero_flag_in){

         // Store descriptor parameters from the snapparam file
         // twojmax sets the highest hyperspherical harmonic angular channel therefore the length/resolution of the bispectrum B
         rfac0 = rfac0_in; 
         twojmax = twojmax_in;
         rmin0 = rmin0_in;
         switch_flag = switch_flag_in;
         bzero_flag = bzero_flag_in;

         // clear neighbour arrays for each central atom
         nmax = 0;
         rij.clear();
         inside.clear();
         wj.clear();
         rcutij.clear();
         ulist_r_ij.clear();
         ulist_i_ij.clear();

         // calculate how many bispectrum coefficients exist for the current twojmax
         compute_ncoeff();
         build_indexlist();

         // allocate the arrays used by the SNA sums. 
         const int jdimpq = twojmax + 2;
         rootpqarray.assign(jdimpq * jdimpq, 0.0);
         cglist.assign(idxcg_max, 0.0);
         ulisttot_r.assign(idxu_max, 0.0);
         ulisttot_i.assign(idxu_max, 0.0);
         zlist_r.assign(idxz_max, 0.0);
         zlist_i.assign(idxz_max, 0.0);
         blist.assign(idxb_max, 0.0);
         dulist_r.assign(idxu_max * 3, 0.0);
         dulist_i.assign(idxu_max * 3, 0.0);
         ylist_r.assign(idxu_max, 0.0);
         ylist_i.assign(idxu_max, 0.0);
         bzero.assign(twojmax + 1, 0.0);

         // bzeroflag subtracts the bispectrum value for an isolated central atom
         if(bzero_flag){
            for(int j = 0; j <= twojmax; j++) bzero[j] = double(j + 1);
         }

         // precompute the coupling coefficients and recursion factors used when building U and coupling U into rotational invariants.
         initialise_clebsch_gordan();
         initialise_rootpqarray();

      }

      // return the index of the block for j1,j2,j
      int snap_sna_t::block_index(const int j1, const int j2, const int j) const{

         const int jdim = twojmax + 1;
         return (j1 * jdim + j2) * jdim + j;

      }

      // return the index of the ulisttot arrays for jj and jju
      int snap_sna_t::u_index(const int jj, const int jju) const{

         return jj * idxu_max + jju;

      }

      // return the index of the dulist arrays for idx and component
      int snap_sna_t::vec_index(const int idx, const int component) const{

         return 3 * idx + component;

      }

      // return the factorial of n as a double
      double snap_sna_t::factorial(const int n) const{

         double value = 1.0;
         for(int i = 2; i <= n; i++) value *= double(i);
         return value;

      }

      // compute the number of bispectrum components for the current twojmax
      void snap_sna_t::compute_ncoeff(){

         int ncount = 0;

         for(int j1 = 0; j1 <= twojmax; j1++){
            for(int j2 = 0; j2 <= j1; j2++){
               // rule: |j1-j2| <= j <= j1+j2 and parity must match (+=2)
               for(int j = j1 - j2; j <= std::min(twojmax, j1 + j2); j += 2){
                  // store symmetry unique bispectrum components
                  if(j >= j1) ncount++;
               }
            }
         }

         ncoeff = ncount;

      }

      // build the packed index tables used by the SNAP sums
      void snap_sna_t::build_indexlist(){

         const int jdim = twojmax + 1;
         idxcg_block.assign(jdim * jdim * jdim, -1);
         idxu_block.assign(jdim, -1);
         idxz_block.assign(jdim * jdim * jdim, -1);
         idxb_block.assign(jdim * jdim * jdim, -1);

         // count and index all Clebsch-Gordan coefficients needed for every valid j1,j2,j block and every m1,m2 pair inside that block
         int idxcg_count = 0;
         for(int j1 = 0; j1 <= twojmax; j1++){
            for(int j2 = 0; j2 <= j1; j2++){
               for(int j = j1 - j2; j <= std::min(twojmax, j1 + j2); j += 2){
                  idxcg_block[block_index(j1,j2,j)] = idxcg_count;
                  for(int m1 = 0; m1 <= j1; m1++){
                     for(int m2 = 0; m2 <= j2; m2++) idxcg_count++;
                  }
               }
            }
         }
         idxcg_max = idxcg_count;

         // U_j has (j+1)^2 complex matrix elements
         int idxu_count = 0;
         for(int j = 0; j <= twojmax; j++){
            idxu_block[j] = idxu_count; // first packed location for channel j
            for(int mb = 0; mb <= j; mb++){
               for(int ma = 0; ma <= j; ma++) idxu_count++;
            }
         }
         idxu_max = idxu_count;

         // B contains just the reduced set of symmetry unique bispectrum triplets
         idxb.clear();
         for(int j1 = 0; j1 <= twojmax; j1++){
            for(int j2 = 0; j2 <= j1; j2++){
               for(int j = j1 - j2; j <= std::min(twojmax, j1 + j2); j += 2){
                  if(j >= j1){
                     b_index_t bindex;
                     bindex.j1 = j1;
                     bindex.j2 = j2;
                     bindex.j = j;
                     idxb_block[block_index(j1,j2,j)] = int(idxb.size()); // map a (j1,j2,j) triple to its B index.
                     idxb.push_back(bindex);
                  }
               }
            }
         }
         idxb_max = int(idxb.size());

         // Z is an intermediate tensor used to form B and Y. 
         // For each angular triple precompute the valid ma/mb ranges so compute_zi() goes through the non-zero coupling terms
         idxz.clear();
         for(int j1 = 0; j1 <= twojmax; j1++){
            for(int j2 = 0; j2 <= j1; j2++){
               for(int j = j1 - j2; j <= std::min(twojmax, j1 + j2); j += 2){
                  idxz_block[block_index(j1,j2,j)] = int(idxz.size());
                  for(int mb = 0; 2 * mb <= j; mb++){
                     for(int ma = 0; ma <= j; ma++){
                        z_index_t zindex;
                        zindex.j1 = j1;
                        zindex.j2 = j2;
                        zindex.j = j;
                        // Clebsch-Gordan selection rules tell us which combinations are non zero
                        zindex.ma1min = std::max(0, (2 * ma - j - j2 + j1) / 2);
                        zindex.ma2max = (2 * ma - j - (2 * zindex.ma1min - j1) + j2) / 2;
                        zindex.na = std::min(j1, (2 * ma - j + j2 + j1) / 2) - zindex.ma1min + 1;
                        zindex.mb1min = std::max(0, (2 * mb - j - j2 + j1) / 2);
                        zindex.mb2max = (2 * mb - j - (2 * zindex.mb1min - j1) + j2) / 2;
                        zindex.nb = std::min(j1, (2 * mb - j + j2 + j1) / 2) - zindex.mb1min + 1;
                        zindex.jju = idxu_block[j] + (j + 1) * mb + ma;
                        idxz.push_back(zindex);
                     }
                  }
               }
            }
         }
         idxz_max = int(idxz.size());

      }

      // Clebsch-Gordan coefficients couple two hyperspherical harmonic expansions into a third angular channel so that the final B components are rotationally invariant
      void snap_sna_t::initialise_clebsch_gordan(){

         int idxcg_count = 0;
         for(int j1 = 0; j1 <= twojmax; j1++){
            for(int j2 = 0; j2 <= j1; j2++){
               for(int j = j1 - j2; j <= std::min(twojmax, j1 + j2); j += 2){
                  for(int m1 = 0; m1 <= j1; m1++){
                     const int aa2 = 2 * m1 - j1;

                     for(int m2 = 0; m2 <= j2; m2++){
                        const int bb2 = 2 * m2 - j2;
                        const int m = (aa2 + bb2 + j) / 2;

                        if(m < 0 || m > j){
                           cglist[idxcg_count] = 0.0;
                           idxcg_count++;
                           continue;
                        }

                        // calculate the Clebsch-Gordan sum for this quantum number combination
                        double sum = 0.0;
                        const int zmin = std::max(0, std::max(-(j - j2 + aa2) / 2, -(j - j1 - bb2) / 2));
                        const int zmax = std::min((j1 + j2 - j) / 2, std::min((j1 - aa2) / 2, (j2 + bb2) / 2));

                        for(int z = zmin; z <= zmax; z++){
                           const int ifac = z % 2 ? -1 : 1;
                           sum += double(ifac) /
                              ( factorial(z) *
                                factorial((j1 + j2 - j) / 2 - z) *
                                factorial((j1 - aa2) / 2 - z) *
                                factorial((j2 + bb2) / 2 - z) *
                                factorial((j - j2 + aa2) / 2 + z) *
                                factorial((j - j1 - bb2) / 2 + z) );
                        }

                        const int cc2 = 2 * m - j;
                        // deltacg and sfaccg are normalisation factors so that the basis is orthonormal when coupling the two U expansions
                        const double dcg = deltacg(j1, j2, j);
                        const double sfaccg = std::sqrt(factorial((j1 + aa2) / 2) *
                                                        factorial((j1 - aa2) / 2) *
                                                        factorial((j2 + bb2) / 2) *
                                                        factorial((j2 - bb2) / 2) *
                                                        factorial((j  + cc2) / 2) *
                                                        factorial((j  - cc2) / 2) *
                                                        double(j + 1));

                        cglist[idxcg_count] = sum * dcg * sfaccg;
                        idxcg_count++;
                     }
                  }
               }
            }
         }

      }

      // return the CG coefficient
      double snap_sna_t::deltacg(const int j1, const int j2, const int j){

         const double sfaccg = factorial((j1 + j2 + j) / 2 + 1);
         return std::sqrt(factorial((j1 + j2 - j) / 2) *
                          factorial((j1 - j2 + j) / 2) *
                          factorial((-j1 + j2 + j) / 2) / sfaccg);

      }

      // initialise repeated sqrt factors in the sums
      void snap_sna_t::initialise_rootpqarray(){

         const int jdimpq = twojmax + 2;
         for(int p = 1; p <= twojmax; p++){
            for(int q = 1; q <= twojmax; q++){
               rootpqarray[p * jdimpq + q] = std::sqrt(double(p) / double(q));
            }
         }

      }

      // increase the neighbour list size and reassign the arrays used to store neighbour info
      void snap_sna_t::grow_rij(const int newnmax){

         if(newnmax <= nmax) return;

         // nmax is the maximum number of neighbours. the current atom can have fewer neighbours 
         // compute_ui() loops through the actual number of neighbours inside the cutoff (j_num)
         nmax = newnmax;
         rij.assign(3 * nmax, 0.0);
         inside.assign(nmax, 0);
         wj.assign(nmax, 0.0);
         rcutij.assign(nmax, 0.0);
         ulist_r_ij.assign(nmax * idxu_max, 0.0);
         ulist_i_ij.assign(nmax * idxu_max, 0.0);

      }

      // zero the total expansion for one central atom before adding the contributions from each neighbour 
      void snap_sna_t::zero_uarraytot(){

         std::fill(ulisttot_r.begin(), ulisttot_r.end(), 0.0);
         std::fill(ulisttot_i.begin(), ulisttot_i.end(), 0.0);

         // add the central atom's self contribution
         for(int j = 0; j <= twojmax; j++){
            int jju = idxu_block[j];
            for(int mb = 0; mb <= j; mb++){
               for(int ma = 0; ma <= j; ma++){
                  if(ma == mb) ulisttot_r[jju] = 1.0;
                  jju++;
               }
            }
         }

      }

      // construct the total wigner expansion for one central atom by adding all neighbours inside their cutoff
      void snap_sna_t::compute_ui(const int jnum){

         zero_uarraytot();

         for(int jj = 0; jj < jnum; jj++){
            const double x = rij[vec_index(jj,0)];
            const double y = rij[vec_index(jj,1)];
            const double z = rij[vec_index(jj,2)];
            const double rsq = x*x + y*y + z*z;
            const double r = std::sqrt(rsq);

            // the 3D neighbour vector is mapped onto a 3-sphere by adding fourth coordinate z0. theta0 controls there the neighbour sits along z0 and is scaled by rfac0 to control the angular resolution of the hypershperical harmonics
            const double theta0 = (r - rmin0) * rfac0 * pi / (rcutij[jj] - rmin0);
            const double z0 = r / std::tan(theta0);

            // compute this neighbour's U matrix
            compute_uarray(x, y, z, z0, r, jj);
            // multiply by cutoff and element weight, then add to the central atom's total neighbour density expansion.
            add_uarraytot(r, jj);
         }

      }

      // compute wigner matrix elements for one neighbour
      void snap_sna_t::compute_uarray(const double x, const double y, const double z,
                                      const double z0, const double r, const int jj){

         const double r0inv = 1.0 / std::sqrt(r*r + z0*z0);
         // a and b are two complex parameters used to represent the point (x,y,z,z0) on the 3-sphere
         // we then build all U^j_{ma,mb}(a,b).
         const double a_r = r0inv * z0;
         const double a_i = -r0inv * z;
         const double b_r = r0inv * y;
         const double b_i = -r0inv * x;

         // j=0 has one basis function, equal to 1.
         ulist_r_ij[u_index(jj,0)] = 1.0;
         ulist_i_ij[u_index(jj,0)] = 0.0;

         const int jdimpq = twojmax + 2;
         for(int j = 1; j <= twojmax; j++){
            int jju = idxu_block[j];
            int jjup = idxu_block[j-1];

            for(int mb = 0; 2 * mb <= j; mb++){
               ulist_r_ij[u_index(jj,jju)] = 0.0;
               ulist_i_ij[u_index(jj,jju)] = 0.0;

               for(int ma = 0; ma < j; ma++){
                  // recursively generate the upper half of the Wigner matrix from the previous angular channel j-1
                  double rootpq = rootpqarray[(j - ma) * jdimpq + (j - mb)];
                  ulist_r_ij[u_index(jj,jju)] +=
                     rootpq * (a_r * ulist_r_ij[u_index(jj,jjup)] +
                               a_i * ulist_i_ij[u_index(jj,jjup)]);
                  ulist_i_ij[u_index(jj,jju)] +=
                     rootpq * (a_r * ulist_i_ij[u_index(jj,jjup)] -
                               a_i * ulist_r_ij[u_index(jj,jjup)]);

                  rootpq = rootpqarray[(ma + 1) * jdimpq + (j - mb)];
                  ulist_r_ij[u_index(jj,jju+1)] =
                     -rootpq * (b_r * ulist_r_ij[u_index(jj,jjup)] +
                                b_i * ulist_i_ij[u_index(jj,jjup)]);
                  ulist_i_ij[u_index(jj,jju+1)] =
                     -rootpq * (b_r * ulist_i_ij[u_index(jj,jjup)] -
                                b_i * ulist_r_ij[u_index(jj,jjup)]);
                  jju++;
                  jjup++;
               }
               jju++;
            }

            // fill the lower half of the Wigner matrix using symmetry
            jju = idxu_block[j];
            jjup = jju + (j + 1) * (j + 1) - 1;
            int mbpar = 1;
            for(int mb = 0; 2 * mb <= j; mb++){
               int mapar = mbpar;
               for(int ma = 0; ma <= j; ma++){
                  if(mapar == 1){
                     ulist_r_ij[u_index(jj,jjup)] = ulist_r_ij[u_index(jj,jju)];
                     ulist_i_ij[u_index(jj,jjup)] = -ulist_i_ij[u_index(jj,jju)];
                  }
                  else{
                     ulist_r_ij[u_index(jj,jjup)] = -ulist_r_ij[u_index(jj,jju)];
                     ulist_i_ij[u_index(jj,jjup)] = ulist_i_ij[u_index(jj,jju)];
                  }
                  mapar = -mapar;
                  jju++;
                  jjup--;
               }
               mbpar = -mbpar;
            }
         }

      }

      // apply the switching function and weight to the neighbour's wigner terms, and accumulate them into the central atom's total wigner expansion
      void snap_sna_t::add_uarraytot(const double r, const int jj){

         // sfac is the smooth radial cutoff f_c(r)
         // wj is the element weight
         // they control how strongly this neighbour contributes to the local density
         const double sfac = compute_sfac(r, rcutij[jj]) * wj[jj];

         for(int j = 0; j <= twojmax; j++){
            int jju = idxu_block[j];
            for(int mb = 0; mb <= j; mb++){
               for(int ma = 0; ma <= j; ma++){
                  // ulisttot is the central atom's full neighbour density
                  // U_total = U_self + sum_neighbours weight*cutoff*U_neighbour
                  ulisttot_r[jju] += sfac * ulist_r_ij[u_index(jj,jju)];
                  ulisttot_i[jju] += sfac * ulist_i_ij[u_index(jj,jju)];
                  jju++;
               }
            }
         }

      }

      // compute intermediate Z terms from the total wigner expansion and the Clebsch-Gordan coefficients
      void snap_sna_t::compute_zi(){

         for(int jjz = 0; jjz < idxz_max; jjz++){
            const int j1 = idxz[jjz].j1;
            const int j2 = idxz[jjz].j2;
            const int ma1min = idxz[jjz].ma1min;
            const int ma2max = idxz[jjz].ma2max;
            const int na = idxz[jjz].na;
            const int mb1min = idxz[jjz].mb1min;
            const int mb2max = idxz[jjz].mb2max;
            const int nb = idxz[jjz].nb;

            const int cgbase = idxcg_block[block_index(j1,j2,idxz[jjz].j)];

            zlist_r[jjz] = 0.0;
            zlist_i[jjz] = 0.0;

            int jju1 = idxu_block[j1] + (j1 + 1) * mb1min;
            int jju2 = idxu_block[j2] + (j2 + 1) * mb2max;
            int icgb = mb1min * (j2 + 1) + mb2max;

            // Z is a Clebsch-Gordan weighted product of two U matrices
            // Calculate U(j1)*U(j2), then loop thorugh and accumulate selection rule allowed ma/mb combos
            for(int ib = 0; ib < nb; ib++){
               double suma1_r = 0.0;
               double suma1_i = 0.0;

               int ma1 = ma1min;
               int ma2 = ma2max;
               int icga = ma1min * (j2 + 1) + ma2max;

               for(int ia = 0; ia < na; ia++){
                  // complex product (u1_r+i u1_i)*(u2_r+i u2_i) is multiplied by the CG coeff for this ma pair
                  suma1_r += cglist[cgbase + icga] *
                     (ulisttot_r[jju1 + ma1] * ulisttot_r[jju2 + ma2] -
                      ulisttot_i[jju1 + ma1] * ulisttot_i[jju2 + ma2]);
                  suma1_i += cglist[cgbase + icga] *
                     (ulisttot_r[jju1 + ma1] * ulisttot_i[jju2 + ma2] +
                      ulisttot_i[jju1 + ma1] * ulisttot_r[jju2 + ma2]);
                  ma1++;
                  ma2--;
                  icga += j2;
               }

               zlist_r[jjz] += cglist[cgbase + icgb] * suma1_r;
               zlist_i[jjz] += cglist[cgbase + icgb] * suma1_i;

               jju1 += j1 + 1;
               jju2 -= j2 + 1;
               icgb += j2;
            }
         }

      }

      // convert the U and Z terms into the linear SNAP bispectrum vector B_i
      void snap_sna_t::compute_bi(){

         for(int jjb = 0; jjb < idxb_max; jjb++){
            const int j1 = idxb[jjb].j1;
            const int j2 = idxb[jjb].j2;
            const int j = idxb[jjb].j;

            int jjz = idxz_block[block_index(j1,j2,j)];
            int jju = idxu_block[j];
            double sumzu = 0.0;

            // B is a contraction of U and Z.
            // Z represents a coupled U(j1)U(j2) product
            // therefore B is a contraction giving the rotationally invariant triple product U(j)*U(j1)*U(j2)
            for(int mb = 0; 2 * mb < j; mb++){
               for(int ma = 0; ma <= j; ma++){
                  sumzu += ulisttot_r[jju] * zlist_r[jjz] +
                           ulisttot_i[jju] * zlist_i[jjz];
                  jjz++;
                  jju++;
               }
            }

            if(j % 2 == 0){
               const int mb = j / 2;
               for(int ma = 0; ma < mb; ma++){
                  sumzu += ulisttot_r[jju] * zlist_r[jjz] +
                           ulisttot_i[jju] * zlist_i[jjz];
                  jjz++;
                  jju++;
               }

               sumzu += 0.5 * (ulisttot_r[jju] * zlist_r[jjz] +
                               ulisttot_i[jju] * zlist_i[jjz]);
            }

            blist[jjb] = 2.0 * sumzu;

            if(bzero_flag) blist[jjb] -= bzero[j]; // subtract the isolated atom density if required
         }

      }

      // compute the linear combination of the bispectrum components to form the final SNAP descriptor vector Y_i for one central atom
      void snap_sna_t::compute_yi(const double* beta){

         // ylist is dE/dU built from the B(U,Z(U)) contraction, once ylist is known, forces only need dU/dR for each neighbour
         std::fill(ylist_r.begin(), ylist_r.end(), 0.0);
         std::fill(ylist_i.begin(), ylist_i.end(), 0.0);

         for(int jjz = 0; jjz < idxz_max; jjz++){
            const int j1 = idxz[jjz].j1;
            const int j2 = idxz[jjz].j2;
            const int j = idxz[jjz].j;
            const int ma1min = idxz[jjz].ma1min;
            const int ma2max = idxz[jjz].ma2max;
            const int na = idxz[jjz].na;
            const int mb1min = idxz[jjz].mb1min;
            const int mb2max = idxz[jjz].mb2max;
            const int nb = idxz[jjz].nb;

            const int cgbase = idxcg_block[block_index(j1,j2,j)];
            double ztmp_r = 0.0;
            double ztmp_i = 0.0;

            int jju1 = idxu_block[j1] + (j1 + 1) * mb1min;
            int jju2 = idxu_block[j2] + (j2 + 1) * mb2max;
            int icgb = mb1min * (j2 + 1) + mb2max;

            // compute a Z-like coupled product (like compute_zi()) and accumulate it into the derivative dE/dU
            for(int ib = 0; ib < nb; ib++){
               double suma1_r = 0.0;
               double suma1_i = 0.0;

               int ma1 = ma1min;
               int ma2 = ma2max;
               int icga = ma1min * (j2 + 1) + ma2max;

               for(int ia = 0; ia < na; ia++){
                  suma1_r += cglist[cgbase + icga] *
                     (ulisttot_r[jju1 + ma1] * ulisttot_r[jju2 + ma2] -
                      ulisttot_i[jju1 + ma1] * ulisttot_i[jju2 + ma2]);
                  suma1_i += cglist[cgbase + icga] *
                     (ulisttot_r[jju1 + ma1] * ulisttot_i[jju2 + ma2] +
                      ulisttot_i[jju1 + ma1] * ulisttot_r[jju2 + ma2]);
                  ma1++;
                  ma2--;
                  icga += j2;
               }

               ztmp_r += cglist[cgbase + icgb] * suma1_r;
               ztmp_i += cglist[cgbase + icgb] * suma1_i;

               jju1 += j1 + 1;
               jju2 -= j2 + 1;
               icgb += j2;
            }

            double betaj = 0.0;
            // the bispectrum only stores non symmetry-equivalent components but the derivative dE/dU must include all permutations of (j1,j2,j), calculate them
            if(j >= j1){
               const int jjb = idxb_block[block_index(j1,j2,j)];
               if(j1 == j){
                  if(j2 == j) betaj = 3.0 * beta[jjb];
                  else betaj = 2.0 * beta[jjb];
               }
               else betaj = beta[jjb];
            }
            else if(j >= j2){
               const int jjb = idxb_block[block_index(j,j2,j1)];
               if(j2 == j) betaj = 2.0 * beta[jjb];
               else betaj = beta[jjb];
            }
            else{
               const int jjb = idxb_block[block_index(j2,j,j1)];
               betaj = beta[jjb];
            }

            if(j1 > j) betaj *= (j1 + 1) / (j + 1.0);

            const int jju = idxz[jjz].jju;
            // accumulate beta weighted Z into ylist=dE/dU for the U matrix element
            ylist_r[jju] += betaj * ztmp_r;
            ylist_i[jju] += betaj * ztmp_i;
         }

      }

      // compute the derivative of the wigner expansion for one neighbour 
      void snap_sna_t::compute_duidrj(const int jj){

         const double x = rij[vec_index(jj,0)];
         const double y = rij[vec_index(jj,1)];
         const double z = rij[vec_index(jj,2)];
         const double rsq = x*x + y*y + z*z;
         const double r = std::sqrt(rsq);
         const double rcut = rcutij[jj];
         const double rscale0 = rfac0 * pi / (rcut - rmin0);
         const double theta0 = (r - rmin0) * rscale0;
         const double cs = std::cos(theta0);
         const double sn = std::sin(theta0);
         // z0 is the fourth coordinate in the SNAP coordinate system
         // dz0/dr is the derivative of z0 wrt neighbour distance r
         const double z0 = r * cs / sn;
         const double dz0dr = z0 / r - (r * rscale0) * (rsq + z0 * z0) / rsq;

         compute_duarray(x, y, z, z0, r, dz0dr, wj[jj], rcut, jj);

      }

      // analytic derivative of the wigner expansion for one neighbour
      void snap_sna_t::compute_duarray(const double x, const double y, const double z,
                                       const double z0, const double r, const double dz0dr,
                                       const double wj_in, const double rcut, const int jj){

         const double rinv = 1.0 / r;
         // unit vector from the central atom to this neighbour, radial derivatives are projected along it (dz0/dr and dz0/dr)
         const double ux = x * rinv;
         const double uy = y * rinv;
         const double uz = z * rinv;

         const double r0inv = 1.0 / std::sqrt(r*r + z0*z0);
         // complex parameters as in compute_uarray()
         const double a_r = z0 * r0inv;
         const double a_i = -z * r0inv;
         const double b_r = y * r0inv;
         const double b_i = -x * r0inv;

         // derivative of the normalisation 1/sqrt(r^2+z0^2)
         const double dr0invdr = -std::pow(r0inv, 3.0) * (r + z0 * dz0dr);

         double dr0inv[3];
         double dz0[3];
         dr0inv[0] = dr0invdr * ux;
         dr0inv[1] = dr0invdr * uy;
         dr0inv[2] = dr0invdr * uz;
         dz0[0] = dz0dr * ux;
         dz0[1] = dz0dr * uy;
         dz0[2] = dz0dr * uz;

         double da_r[3], da_i[3], db_r[3], db_i[3];
         for(int k = 0; k < 3; k++){
            // derivatives of a and b in cartesian coordinates
            da_r[k] = dz0[k] * r0inv + z0 * dr0inv[k];
            da_i[k] = -z * dr0inv[k];
            db_r[k] = y * dr0inv[k];
            db_i[k] = -x * dr0inv[k];
         }

         // add the derivative of -i*z, -i*x, and y with respect to the corresponding cartesian components
         da_i[2] += -r0inv;
         db_i[0] += -r0inv;
         db_r[1] += r0inv;

         std::fill(dulist_r.begin(), dulist_r.end(), 0.0);
         std::fill(dulist_i.begin(), dulist_i.end(), 0.0);

         const int jdimpq = twojmax + 2;
         for(int j = 1; j <= twojmax; j++){
            int jju = idxu_block[j];
            int jjup = idxu_block[j-1];

            for(int mb = 0; 2 * mb <= j; mb++){
               for(int ma = 0; ma < j; ma++){
                  double rootpq = rootpqarray[(j - ma) * jdimpq + (j - mb)];
                  for(int k = 0; k < 3; k++){
                     // differentiate the recursion relation for U^j_{ma,mb}(a,b) with respect to the cartesian coordinates of the neighbour
                     // Each term has derivative of the basis parameter (da/db) times old U plus the basis parameter times old dU
                     dulist_r[vec_index(jju,k)] +=
                        rootpq * (da_r[k] * ulist_r_ij[u_index(jj,jjup)] +
                                  da_i[k] * ulist_i_ij[u_index(jj,jjup)] +
                                  a_r * dulist_r[vec_index(jjup,k)] +
                                  a_i * dulist_i[vec_index(jjup,k)]);
                     dulist_i[vec_index(jju,k)] +=
                        rootpq * (da_r[k] * ulist_i_ij[u_index(jj,jjup)] -
                                  da_i[k] * ulist_r_ij[u_index(jj,jjup)] +
                                  a_r * dulist_i[vec_index(jjup,k)] -
                                  a_i * dulist_r[vec_index(jjup,k)]);
                  }

                  rootpq = rootpqarray[(ma + 1) * jdimpq + (j - mb)];
                  for(int k = 0; k < 3; k++){
                     dulist_r[vec_index(jju+1,k)] =
                        -rootpq * (db_r[k] * ulist_r_ij[u_index(jj,jjup)] +
                                   db_i[k] * ulist_i_ij[u_index(jj,jjup)] +
                                   b_r * dulist_r[vec_index(jjup,k)] +
                                   b_i * dulist_i[vec_index(jjup,k)]);
                     dulist_i[vec_index(jju+1,k)] =
                        -rootpq * (db_r[k] * ulist_i_ij[u_index(jj,jjup)] -
                                   db_i[k] * ulist_r_ij[u_index(jj,jjup)] +
                                   b_r * dulist_i[vec_index(jjup,k)] -
                                   b_i * dulist_r[vec_index(jjup,k)]);
                  }

                  jju++;
                  jjup++;
               }
               jju++;
            }

            jju = idxu_block[j];
            jjup = jju + (j + 1) * (j + 1) - 1;
            int mbpar = 1;
            // fill the derivative lower half using the same symmetry used for the U matrix itself
            for(int mb = 0; 2 * mb <= j; mb++){
               int mapar = mbpar;
               for(int ma = 0; ma <= j; ma++){
                  if(mapar == 1){
                     for(int k = 0; k < 3; k++){
                        dulist_r[vec_index(jjup,k)] = dulist_r[vec_index(jju,k)];
                        dulist_i[vec_index(jjup,k)] = -dulist_i[vec_index(jju,k)];
                     }
                  }
                  else{
                     for(int k = 0; k < 3; k++){
                        dulist_r[vec_index(jjup,k)] = -dulist_r[vec_index(jju,k)];
                        dulist_i[vec_index(jjup,k)] = dulist_i[vec_index(jju,k)];
                     }
                  }
                  mapar = -mapar;
                  jju++;
                  jjup--;
               }
               mbpar = -mbpar;
            }
         }

         const double sfac = compute_sfac(r, rcut) * wj_in;
         const double dsfac = compute_dsfac(r, rcut) * wj_in;

         for(int j = 0; j <= twojmax; j++){
            int jju = idxu_block[j];
            for(int mb = 0; 2 * mb <= j; mb++){
               for(int ma = 0; ma <= j; ma++){
                  // product rule: d(sfac * U)/dR = dsfac/dr * rhat * U + sfac * dU/dR.
                  dulist_r[vec_index(jju,0)] = dsfac * ulist_r_ij[u_index(jj,jju)] * ux +
                                               sfac * dulist_r[vec_index(jju,0)];
                  dulist_i[vec_index(jju,0)] = dsfac * ulist_i_ij[u_index(jj,jju)] * ux +
                                               sfac * dulist_i[vec_index(jju,0)];
                  dulist_r[vec_index(jju,1)] = dsfac * ulist_r_ij[u_index(jj,jju)] * uy +
                                               sfac * dulist_r[vec_index(jju,1)];
                  dulist_i[vec_index(jju,1)] = dsfac * ulist_i_ij[u_index(jj,jju)] * uy +
                                               sfac * dulist_i[vec_index(jju,1)];
                  dulist_r[vec_index(jju,2)] = dsfac * ulist_r_ij[u_index(jj,jju)] * uz +
                                               sfac * dulist_r[vec_index(jju,2)];
                  dulist_i[vec_index(jju,2)] = dsfac * ulist_i_ij[u_index(jj,jju)] * uz +
                                               sfac * dulist_i[vec_index(jju,2)];
                  jju++;
               }
            }
         }

      }

      // finish the force calculation for one neighbour (dE_i/dR_ij) then add it to the central atom and subtract from the neighbour
      void snap_sna_t::compute_deidrj(double* dedr){

         dedr[0] = 0.0;
         dedr[1] = 0.0;
         dedr[2] = 0.0;

         for(int j = 0; j <= twojmax; j++){
            int jju = idxu_block[j];

            for(int mb = 0; 2 * mb < j; mb++){
               for(int ma = 0; ma <= j; ma++){
                  const double yr = ylist_r[jju];
                  const double yi = ylist_i[jju];
                  for(int k = 0; k < 3; k++){
                     // chain rule: dE/dR_k = sum_U (dU/dR_k) * (dE/dU) - handle Re and Im parts of the complex U matrix separately
                     dedr[k] += dulist_r[vec_index(jju,k)] * yr +
                                dulist_i[vec_index(jju,k)] * yi;
                  }
                  jju++;
               }
            }

            if(j % 2 == 0){
               const int mb = j / 2;
               for(int ma = 0; ma < mb; ma++){
                  const double yr = ylist_r[jju];
                  const double yi = ylist_i[jju];
                  for(int k = 0; k < 3; k++){
                     dedr[k] += dulist_r[vec_index(jju,k)] * yr +
                                dulist_i[vec_index(jju,k)] * yi;
                  }
                  jju++;
               }

               const double yr = ylist_r[jju];
               const double yi = ylist_i[jju];
               for(int k = 0; k < 3; k++){
                  dedr[k] += 0.5 * (dulist_r[vec_index(jju,k)] * yr +
                                    dulist_i[vec_index(jju,k)] * yi);
               }
            }
         }

         dedr[0] *= 2.0;
         dedr[1] *= 2.0;
         dedr[2] *= 2.0;

      }

      // default cosine cutoff for SNAP
      double snap_sna_t::compute_sfac(const double r, const double rcut){

         if(switch_flag == 0) return 1.0;
         if(r <= rmin0) return 1.0;
         if(r > rcut) return 0.0;

         // it smoothly cuts off the neighbour density from 1 to 0 between rmin0 and rcut which removes discontinuities in energy when neighbours enter or leave the descriptor environment
         const double rcutfac = pi / (rcut - rmin0);
         return 0.5 * (std::cos((r - rmin0) * rcutfac) + 1.0);

      }

      // derivative of the cosine cutoff
      double snap_sna_t::compute_dsfac(const double r, const double rcut){

         if(switch_flag == 0) return 0.0;
         if(r <= rmin0) return 0.0;
         if(r > rcut) return 0.0;

         // d f_c(r) / dr used in the force calculation
         const double rcutfac = pi / (rcut - rmin0);
         return -0.5 * std::sin((r - rmin0) * rcutfac) * rcutfac;

      }


      // this is visible to the rest of Vampire and is used to store the SNAP potential parameters and coefficients
      snap_potential_t::snap_potential_t() :
         initialised(false),
         debug_enabled(false),
         debug_force_calls(0),
         debug_max_force_calls(5),
         rcutfac(0.0),
         rfac0(0.99363),
         rmin0(0.0),
         rcutmax(0.0),
         twojmax(0),
         switchflag(1),
         bzeroflag(1),
         ncoeff(0)
      {
      }

      // set the SNAP coefficient file
      void snap_potential_t::set_coeff_filename(const std::string& filename){

         coeff_filename = filename;

      }

      // set the SNAP parameter file
      void snap_potential_t::set_param_filename(const std::string& filename){

         param_filename = filename;

      }

      void snap_potential_t::set_debug(const bool debug){

         debug_enabled = debug;
         debug_force_calls = 0;

      }

      bool snap_potential_t::debug_this_call() const{

         return debug_enabled && debug_force_calls < debug_max_force_calls;

      }

      int snap_potential_t::debug_call_number() const{

         return debug_force_calls + 1;

      }

      void snap_potential_t::increment_debug_force_calls(){

         debug_force_calls++;

      }

      // return the max cutoff for SNAP
      double snap_potential_t::cutoff() const{

         return rcutmax;

      }

      // return the number of coefficients in the SNAP linear model (not including the intercept)
      int snap_potential_t::number_of_coefficients() const{

         return ncoeff;

      }

      // return the number of elements in the SNAP potential (currently only one element is supported)
      int snap_potential_t::number_of_elements() const{

         return int(elements.size());

      }

      std::string snap_potential_t::element_name(const int element) const{

         if(element < 0 || element >= int(elements.size())) return "";
         return elements[element].name;

      }

      // initialise the SNAP potential by reading the coefficient and parameter files, and setting up the data structures for the descriptor calcs
      void snap_potential_t::initialise(const int num_materials){

         if(coeff_filename.empty()){
            err::zexit("SNAP potential selected but spin-lattice:snap-coeff-file was not set");
         }

         if(param_filename.empty()){
            err::zexit("SNAP potential selected but spin-lattice:snap-param-file was not set");
         }

         // Read the SNAP files 
         // The coefficient file defines the fitted linear model E = beta0 + sum_k beta_k B_k
         // The parameter file defines how the B_k descriptors are constructed
         read_files();
         // map Vampire material ids to SNAP element ids
         build_material_map(num_materials); // SNAP only supports one element at the moment so all mapped to the first element

         if(elements.size() == 0 || elements[0].coefficients.size() <= 1){
            err::zexit("Incorrect SNAP coefficient file");
         }
         // ncoeff counts only the descriptor coefficients beta_1...beta_n (not the intercept beta_0)
         // there are ncoeff + 1 coefficients in the file, including the intercept
         ncoeff = int(elements[0].coefficients.size()) - 1;

         // initialise descriptor data structures after reading twojmax so the number and order of B components matches the coefficient file
         sna.initialise(rfac0, twojmax, rmin0, switchflag, bzeroflag);
         if(ncoeff != sna.ncoeff){
            err::zexit("SNAP coefficient count does not match twojmax in SNAP parameter file");
         }

         rcutmax = 0.0;
         for(unsigned int element = 0; element < elements.size(); element++){
            // multi-element SNAP pair cutoff from element radii and rcutfac
            // currently single element in vampire so this becomes rcutmax = 2*radius_Fe*rcutfac.
            rcutmax = std::max(rcutmax, 2.0 * elements[element].radius * rcutfac);
         }

         initialised = true;
         print_info();
         if(debug_enabled){
            std::cout<<"SNAP debug output enabled for the first " <<debug_max_force_calls<<" force calls"<<std::endl;
         }

      }

      // read the SNAP coefficient and parameter files
      void snap_potential_t::read_files(){

         read_coeff_file();
         read_param_file();

      }

      // read the SNAP coefficient file and store the coefficients for each element
      void snap_potential_t::read_coeff_file(){

         std::ifstream file(coeff_filename.c_str());
         if(!file){
            err::zexit("Cannot open SNAP coefficient file \"" + coeff_filename + "\"");
         }

         std::vector<std::string> header = read_words(file);
         if(header.size() != 2){
            err::zexit("Incorrect format in SNAP coefficient file \"" + coeff_filename + "\"");
         }

         // first non-comment line is the number of element blocks and number of coefficients in each block, including the intercept.
         const int nelements = inumeric(header[0], coeff_filename);
         const int ncoeffall = inumeric(header[1], coeff_filename);
         if(nelements <= 0 || ncoeffall <= 0){
            err::zexit("Incorrect header in SNAP coefficient file \"" + coeff_filename + "\"");
         }

         // read the element blocks and store the coefficients for each element
         elements.clear();
         elements.reserve(nelements);

         for(int element = 0; element < nelements; element++){
            std::vector<std::string> element_header = read_words(file);
            if(element_header.size() != 3){
               err::zexit("Incorrect element block in SNAP coefficient file \"" + coeff_filename + "\"");
            }

            snap_element_t snap_element;
            snap_element.name = element_header[0];
            // radius and weight are descriptor parameters used by SNAP to form neighbour density, not to be confused with material mass or the physical radius
            snap_element.radius = numeric(element_header[1], coeff_filename);
            snap_element.weight = numeric(element_header[2], coeff_filename);
            snap_element.coefficients.resize(ncoeffall, 0.0);

            for(int icoeff = 0; icoeff < ncoeffall; icoeff++){
               std::vector<std::string> coefficient_line = read_words(file);
               if(coefficient_line.size() != 1){
                  err::zexit("Incorrect coefficient block in SNAP coefficient file \"" + coeff_filename + "\"");
               }
               // coefficient[0] is beta0 which is the per-atom energy intercept
               // coefficient[k+1] multiplies bispectrum component B_k
               snap_element.coefficients[icoeff] = numeric(coefficient_line[0], coeff_filename);
            }

            elements.push_back(snap_element);
         }

      }

      void snap_potential_t::read_param_file(){

         std::ifstream file(param_filename.c_str());
         if(!file){
            err::zexit("Cannot open SNAP parameter file \"" + param_filename + "\"");
         }

         // default settings for the SNAP parameter file
         rcutfac = 0.0;
         rfac0 = 0.99363;
         rmin0 = 0.0;
         twojmax = -1;
         switchflag = 1;
         bzeroflag = 1;

         int chemflag = 0;
         int bnormflag = 0;
         int wselfallflag = 0;
         int switchinnerflag = 0;

         std::string line;
         while(std::getline(file, line)){
            std::vector<std::string> words = words_from_line(line);
            if(words.size() == 0) continue;
            if(words.size() != 2){
               err::zexit("Incorrect format in SNAP parameter file \"" + param_filename + "\"");
            }

            const std::string key = words[0];
            const std::string value = words[1];

            // current implementation supports linear, single element SNAP so reject any unsupported settings in the parameter file 
            if(key == "rcutfac") rcutfac = numeric(value, param_filename);
            else if(key == "twojmax") twojmax = inumeric(value, param_filename);
            else if(key == "rfac0") rfac0 = numeric(value, param_filename);
            else if(key == "rmin0") rmin0 = numeric(value, param_filename);
            else if(key == "switchflag") switchflag = inumeric(value, param_filename);
            else if(key == "bzeroflag") bzeroflag = inumeric(value, param_filename);
            else if(key == "quadraticflag"){
               if(inumeric(value, param_filename) != 0){
                  err::zexit("SNAP quadraticflag is not supported");
               }
            }
            else if(key == "chemflag") chemflag = inumeric(value, param_filename);
            else if(key == "bnormflag") bnormflag = inumeric(value, param_filename);
            else if(key == "wselfallflag") wselfallflag = inumeric(value, param_filename);
            else if(key == "switchinnerflag") switchinnerflag = inumeric(value, param_filename);
            else if(key == "chunksize") continue;
            else if(key == "parallelthresh") continue;
            else{
               err::zexit("Unknown keyword \"" + key + "\" in SNAP parameter file \"" + param_filename + "\"");
            }
         }

         if(rcutfac <= 0.0 || twojmax < 0){
            err::zexit("SNAP parameter file must set rcutfac and twojmax");
         }

         if(chemflag != 0){
            err::zexit("SNAP chemflag is not supported");
         }

         if(bnormflag != 0){
            err::zexit("SNAP bnormflag is not supported");
         }

         if(wselfallflag != 0){
            err::zexit("SNAP wselfallflag is not supported");
         }

         if(switchinnerflag != 0){
            err::zexit("SNAP switchinnerflag is not supported");
         }

      }

      void snap_potential_t::build_material_map(const int num_materials){

         material_to_element.assign(num_materials, 0);

         // currently supports a single element so all Vampire mat ids are mapped to that element
         if(elements.size() != 1){
            err::zexit("SNAP potential currently supports single-element SNAP coefficient files only");
         }

      }

      // print the SNAP potential information to the terminal
      void snap_potential_t::print_info() const{

         std::cout<<"SNAP potential information:"<<std::endl;
         std::cout<<"  Coefficient file: "<<coeff_filename<<std::endl;
         std::cout<<"  Parameter file: "<<param_filename<<std::endl;
         std::cout<<"  rcutfac: "<<rcutfac<<std::endl;
         std::cout<<"  rfac0: "<<rfac0<<std::endl;
         std::cout<<"  rmin0: "<<rmin0<<std::endl;
         std::cout<<"  twojmax: "<<twojmax<<std::endl;
         std::cout<<"  switchflag: "<<switchflag<<std::endl;
         std::cout<<"  bzeroflag: "<<bzeroflag<<std::endl;
         std::cout<<"  cutoff: "<<rcutmax<<std::endl;
         std::cout<<"  descriptor coefficients: "<<ncoeff<<std::endl;

         for(unsigned int element = 0; element < elements.size(); element++){
            std::cout<<"  Element "<<element<<": "<<elements[element].name<<std::endl;
            std::cout<<"    radius: "<<elements[element].radius<<std::endl;
            std::cout<<"    weight: "<<elements[element].weight<<std::endl;
            std::cout<<"    coefficients including intercept:"<<std::endl;
            for(unsigned int coefficient = 0; coefficient < elements[element].coefficients.size(); coefficient++){
               std::cout<<"      "<<coefficient<<" "
                        <<std::setprecision(17)
                        <<elements[element].coefficients[coefficient]
                        <<std::setprecision(6)<<std::endl;
            }
         }

      }

      // construct the short neighbour list for one atom, and store the neighbour distances and weights in the SNAP data structure
      int snap_potential_t::build_short_neighbour_list(const int atom,
                                                       const std::vector<int>& neighbour_list_start_index,
                                                       const std::vector<int>& neighbour_list_end_index,
                                                       const std::vector<int>& type_array,
                                                       const std::vector<int>& neighbour_list_array,
                                                       const std::vector<double>& x_coord_array,
                                                       const std::vector<double>& y_coord_array,
                                                       const std::vector<double>& z_coord_array){

         const int imat = type_array[atom];
         if(imat < 0 || imat >= int(material_to_element.size())){
            err::zexit("SNAP material index is outside the element map");
         }

         const int ielem = material_to_element[imat];
         const double radi = elements[ielem].radius;
         const double xtmp = x_coord_array[atom];
         const double ytmp = y_coord_array[atom];
         const double ztmp = z_coord_array[atom];

         // start from Vampire's spin-lattice neighbour list, then filter it based on the true SNAP distance cutoff
         // this means Vampire's neighbour list must be sufficiently large enough that it contains every atom inside rcutij
         const int nbr_start = neighbour_list_start_index[atom];
         const int nbr_end = neighbour_list_end_index[atom] + 1;
         sna.grow_rij(nbr_end - nbr_start);

         int ninside = 0;
         for(int jj = nbr_start; jj < nbr_end; jj++){
            const int j = neighbour_list_array[jj];
            if(j == atom) continue;

            const int jmat = type_array[j];
            if(jmat < 0 || jmat >= int(material_to_element.size())){
               err::zexit("SNAP neighbour material index is outside the element map");
            }

            const int jelem = material_to_element[jmat];
            // vector from central atom to the neighbour
            double delx = x_coord_array[j] - xtmp;
            double dely = y_coord_array[j] - ytmp;
            double delz = z_coord_array[j] - ztmp;

            // PBC wrapping chooses the nearest image before the cutoff test
            delx = sld::PBC_wrap(delx, cs::system_dimensions[0], cs::pbc[0]);
            dely = sld::PBC_wrap(dely, cs::system_dimensions[1], cs::pbc[1]);
            delz = sld::PBC_wrap(delz, cs::system_dimensions[2], cs::pbc[2]);

            const double rsq = delx*delx + dely*dely + delz*delz;
            // element-pair SNAP cutoff, for single-element this is the same for every pair
            const double rcutij = (radi + elements[jelem].radius) * rcutfac;

            if(rsq < rcutij * rcutij && rsq > 1.0e-20){
               // store the neighbour environment fed into the SNAP descriptor calculation
               // the neighbour list now only includes neighbours inside the cutoff
               // the order matches vampires neighbour list and defines jj indices used when differentiating neighbour contributions
               sna.rij[3 * ninside + 0] = delx;
               sna.rij[3 * ninside + 1] = dely;
               sna.rij[3 * ninside + 2] = delz;
               sna.inside[ninside] = j;
               sna.wj[ninside] = elements[jelem].weight;
               sna.rcutij[ninside] = rcutij;
               ninside++;
            }
         }

         return ninside;

      }

      // compute the bispectrum for a range of atoms, using the short neighbour list and the SNAP data structure then store it in the output vector
      void snap_potential_t::compute_bispectrum(const int start_index,
                                                const int end_index,
                                                const std::vector<int>& neighbour_list_start_index,
                                                const std::vector<int>& neighbour_list_end_index,
                                                const std::vector<int>& type_array,
                                                const std::vector<int>& neighbour_list_array,
                                                const std::vector<double>& x_coord_array,
                                                const std::vector<double>& y_coord_array,
                                                const std::vector<double>& z_coord_array,
                                                std::vector<double>& bispectrum){

         for(int i = start_index; i < end_index; i++){
            const int local_atom = i - start_index;
            // build the central atom environment first, ninside is the number of neighbours that actually lie inside the element-pair cutoff.
            const int ninside = build_short_neighbour_list(i,
                                                           neighbour_list_start_index,
                                                           neighbour_list_end_index,
                                                           type_array,
                                                           neighbour_list_array,
                                                           x_coord_array,
                                                           y_coord_array,
                                                           z_coord_array);

            // Ui is the neighbour density expansion, Zi is the Clebsch-Gordan coupled intermediate and Bi is the rotationally invariant bispectrum descriptor vector for this central atom
            sna.compute_ui(ninside);
            sna.compute_zi();
            sna.compute_bi();

            for(int icoeff = 0; icoeff < ncoeff; icoeff++){
               // store B_i,k in a flat array so later energy and debug code can index using local_atom*ncoeff + icoeff
               bispectrum[local_atom * ncoeff + icoeff] = sna.blist[icoeff];
            }
         }

      }

      // compute the derivative of the energy wrt the bispectrum for atoms
      void snap_potential_t::compute_beta(const int start_index,
                                          const int end_index,
                                          const std::vector<int>& type_array,
                                          std::vector<double>& beta) const{

         for(int i = start_index; i < end_index; i++){
            const int local_atom = i - start_index;
            const int imat = type_array[i];
            if(imat < 0 || imat >= int(material_to_element.size())){
               err::zexit("SNAP material index is outside the element map");
            }
            const int ielem = material_to_element[imat];
            const std::vector<double>& coeffi = elements[ielem].coefficients;

            // linear SNAP energy is E_i = beta0 + sum_k beta_k B_i,k, so dE_i/dB_i,k is the fitted coefficient beta_k
            for(int icoeff = 0; icoeff < ncoeff; icoeff++){
               beta[local_atom * ncoeff + icoeff] = coeffi[icoeff + 1];
            }
         }

      }

      // compute the energy for one atom
      double snap_potential_t::compute_energy(const int atom,
                                              const int local_atom,
                                              const std::vector<int>& type_array,
                                              const std::vector<double>& bispectrum) const{

         const int imat = type_array[atom];
         if(imat < 0 || imat >= int(material_to_element.size())){
            err::zexit("SNAP material index is outside the element map");
         }
         const int ielem = material_to_element[imat];
         const std::vector<double>& coeffi = elements[ielem].coefficients;

         // the first coefficient is the constant per-atom intercept
         double energy = coeffi[0];
         for(int icoeff = 0; icoeff < ncoeff; icoeff++){
            // add the linear contribution from each bispectrum descriptor
            energy += coeffi[icoeff + 1] * bispectrum[local_atom * ncoeff + icoeff];
         }

         return energy;

      }

      // print the force debugging for SNAP
      void snap_potential_t::print_debug_summary(const int call,
                                                 const int start_index,
                                                 const int end_index,
                                                 const std::vector<int>& neighbour_counts,
                                                 const std::vector<double>& initial_forces_x,
                                                 const std::vector<double>& initial_forces_y,
                                                 const std::vector<double>& initial_forces_z,
                                                 const std::vector<double>& forces_array_x,
                                                 const std::vector<double>& forces_array_y,
                                                 const std::vector<double>& forces_array_z) const{

         sld::internal::print_force_debug_summary("SNAP",
                                                  call,
                                                  start_index, end_index,
                                                  neighbour_counts,
                                                  initial_forces_x,
                                                  initial_forces_y,
                                                  initial_forces_z,
                                                  forces_array_x,
                                                  forces_array_y,
                                                  forces_array_z);

      }

      // compute the forces for atoms
      void snap_potential_t::compute_forces(const int start_index,
                                            const int end_index,
                                            const std::vector<int>& neighbour_list_start_index,
                                            const std::vector<int>& neighbour_list_end_index,
                                            const std::vector<int>& type_array,
                                            const std::vector<int>& neighbour_list_array,
                                            const std::vector<double>& x_coord_array,
                                            const std::vector<double>& y_coord_array,
                                            const std::vector<double>& z_coord_array,
                                            std::vector<double>& forces_array_x,
                                            std::vector<double>& forces_array_y,
                                            std::vector<double>& forces_array_z,
                                            std::vector<double>& potential_eng,
                                            const bool allow_debug_output){

         const int natoms = end_index - start_index;
         std::vector<double> beta(natoms * ncoeff, 0.0);
         std::vector<double> bispectrum(natoms * ncoeff, 0.0);
         const bool debug_this_call = allow_debug_output && this->debug_this_call();
         std::vector<int> neighbour_counts;
         std::vector<double> initial_forces_x;
         std::vector<double> initial_forces_y;
         std::vector<double> initial_forces_z;

         // for the debugging, save the force arrays for comparison
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

         // first compute B_i,k for every central atom in this range
         // then compute beta_i,k = dE_i/dB_i,k from the fitted coefficients
         // finally go back to each each atom and differentiate E_i with respect to each neighbour coordinate to get forces from that neighbour 
         compute_bispectrum(start_index, end_index,
                            neighbour_list_start_index, neighbour_list_end_index,
                            type_array, neighbour_list_array,
                            x_coord_array, y_coord_array, z_coord_array,
                            bispectrum);
         compute_beta(start_index, end_index, type_array, beta);

         for(int i = start_index; i < end_index; i++){
            const int local_atom = i - start_index;
            // once B is known, calculate the per-atom SNAP energy
            // In snapzbl the per-atom SNAP_energy is stored in potential_eng[i] so that the ZBL potential overlay can add its pair energy to this value
            potential_eng[i] = compute_energy(i, local_atom, type_array, bispectrum);

            // rebuild the same neighbour environment for calculating derivatives
            const int ninside = build_short_neighbour_list(i,
                                                           neighbour_list_start_index,
                                                           neighbour_list_end_index,
                                                           type_array,
                                                           neighbour_list_array,
                                                           x_coord_array,
                                                           y_coord_array,
                                                           z_coord_array);
            if(debug_this_call){
               neighbour_counts[local_atom] = ninside;
            }

            // compute U_i for this central atom again, then form Y_i=dE_i/dU_i by contracting beta with the Z intermediates. 
            // Y_i is calcualted so that each neighbour force can be computed as dU/dR*dE/dU
            sna.compute_ui(ninside);
            sna.compute_yi(&beta[local_atom * ncoeff]);

            // for the neighbours of the central atom within cutoff:
            // compute Fij = dE_i/dR_ij. 
            // since R_ij is neighbour position relative to the central atom, we add this contribution to the central atom and subtract from the neighbour to conserve action/reaction
            for(int jj = 0; jj < ninside; jj++){
               double fij[3];
               const int j = sna.inside[jj];

               // differentiate this neighbour's U contribution with respect to its x/y/z displacement from the central atom
               sna.compute_duidrj(jj);
               // then contract dU/dR with Y=dE/dU to obtain the force vector
               sna.compute_deidrj(fij);

               forces_array_x[i] += fij[0];
               forces_array_y[i] += fij[1];
               forces_array_z[i] += fij[2];
               forces_array_x[j] -= fij[0];
               forces_array_y[j] -= fij[1];
               forces_array_z[j] -= fij[2];
            }
         }

         if(debug_this_call){
            print_debug_summary(debug_force_calls + 1,
                                start_index, end_index,
                                neighbour_counts,
                                initial_forces_x,
                                initial_forces_y,
                                initial_forces_z,
                                forces_array_x,
                                forces_array_y,
                                forces_array_z);
            increment_debug_force_calls();
         }

      }

   } // end of internal namespace

} // end of sld namespace
