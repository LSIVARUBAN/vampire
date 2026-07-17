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

#ifndef SLD_SNAP_H_
#define SLD_SNAP_H_

// C++ standard library headers
#include <string>
#include <vector>

namespace sld{

   namespace internal{

      // stores params read from one element block in the snapcoeff file.
      // radius is used to form the element pair cutoff rcut_ij = (radius_i + radius_j) * rcutfac, and weight to scale that neighbour's contribution to the local density expansion
      struct snap_element_t{
         std::string name;
         double radius;
         double weight;
         std::vector<double> coefficients;
      };

      // Bispectrum descriptor class
      // For one central atom it stores the neighbour vectors, expands them in 4D hyperspherical harmonics U, couples U terms with Clebsch-Gordan coefficients to build Z, then contracts U and Z into the rotationally invariant bispectrum B
      class snap_sna_t{

      public:

         snap_sna_t();

         void initialise(const double rfac0_in,
                         const int twojmax_in,
                         const double rmin0_in,
                         const int switch_flag_in,
                         const int bzero_flag_in);

         void grow_rij(const int newnmax);
         void compute_ui(const int jnum);
         void compute_zi();
         void compute_bi();
         void compute_yi(const double* beta);
         void compute_duidrj(const int jj);
         void compute_deidrj(double* dedr);

         int ncoeff;
         int nmax;
         int twojmax;

         std::vector<double> blist;  // final B_k bispectrum components for the current central atom
         std::vector<double> rij;    // packed neighbour vectors relative to the current central atom
         std::vector<int> inside;    // atom indices of neighbours inside the SNAP cutoff
         std::vector<double> wj;     // element weights for each neighbour contribution
         std::vector<double> rcutij; // element-pair cutoffs for each neighbour

      private:

         // Packed indexing for one Z term. Z couples two U expansions, j1 and j2 into j
         // ma and mb ranges avoid looping over coefficient combinations that Clebsch-Gordan selection rules make zero.
         struct z_index_t{
            int j1;
            int j2;
            int j;
            int ma1min;
            int ma2max;
            int mb1min;
            int mb2max;
            int na;
            int nb;
            int jju;
         };

         // One bispectrum component B_{j1,j2,j}
         struct b_index_t{
            int j1;
            int j2;
            int j;
         };

         double rfac0;
         double rmin0;
         int switch_flag;
         int bzero_flag;

         int idxcg_max;
         int idxu_max;
         int idxz_max;
         int idxb_max;

         std::vector<z_index_t> idxz;
         std::vector<b_index_t> idxb;
         std::vector<int> idxcg_block;
         std::vector<int> idxu_block;
         std::vector<int> idxz_block;
         std::vector<int> idxb_block;

         std::vector<double> rootpqarray; // sqrt factors used in the recursive U construction
         std::vector<double> cglist;      // Clebsch-Gordan coefficients
         std::vector<double> ulisttot_r;  // real part of total neighbour density expansion U
         std::vector<double> ulisttot_i;  // imaginary part of total neighbour density expansion U
         std::vector<double> ulist_r_ij;  // real U contribution from each individual neighbour
         std::vector<double> ulist_i_ij;  // imaginary U contribution from each individual neighbour
         std::vector<double> zlist_r;     // real intermediate Z terms used to form B and Y
         std::vector<double> zlist_i;     // imaginary intermediate Z terms
         std::vector<double> dulist_r;    // real derivative dU/dR for one neighbour
         std::vector<double> dulist_i;    // imaginary derivative dU/dR for one neighbour
         std::vector<double> ylist_r;     // real derivative dE/dU after contracting beta with Z
         std::vector<double> ylist_i;     // imaginary derivative dE/dU after contracting beta with Z
         std::vector<double> bzero;       // isolated atom B offsets subtracted if bzeroflag is enabled

         void build_indexlist();
         void initialise_clebsch_gordan();
         void initialise_rootpqarray();
         void compute_ncoeff();
         void zero_uarraytot();
         void add_uarraytot(const double r, const int jj);
         void compute_uarray(const double x, const double y, const double z,
                             const double z0, const double r, const int jj);
         void compute_duarray(const double x, const double y, const double z,
                              const double z0, const double r, const double dz0dr,
                              const double wj_in, const double rcut, const int jj);
         double deltacg(const int j1, const int j2, const int j);
         double compute_sfac(const double r, const double rcut);
         double compute_dsfac(const double r, const double rcut);
         int block_index(const int j1, const int j2, const int j) const;
         int u_index(const int jj, const int jju) const;
         int vec_index(const int idx, const int component) const;
         double factorial(const int n) const;
      };

      class snap_potential_t{

      public:

         snap_potential_t();

         void initialise(const int num_materials);

         void compute_forces(const int start_index,
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
                             const bool allow_debug_output = true);

         void set_coeff_filename(const std::string& filename);
         void set_param_filename(const std::string& filename);
         void set_debug(const bool debug);
         bool debug_this_call() const;
         int debug_call_number() const;
         void increment_debug_force_calls();

         // SNAP files contain the cutoff information
         double cutoff() const;
         int number_of_coefficients() const;
         int number_of_elements() const;
         std::string element_name(const int element) const;

      private:

         bool initialised;
         bool debug_enabled;
         int debug_force_calls;
         int debug_max_force_calls;
         std::string coeff_filename; // coefficients
         std::string param_filename; // descriptor settings
         std::vector<snap_element_t> elements; // SNAP elements *supports only one currently*
         std::vector<int> material_to_element; // maps Vampire material ids to SNAP element ids
         snap_sna_t sna; // Bispectrum descriptor class

         double rcutfac;  // multiplier used with element radii to form the pair cutoff
         double rfac0;    // angular scaling factor for 3D-4D mapping
         double rmin0;    // inner radius for the SNAP cutoff and coordinate mapping
         double rcutmax;  // maximum element-pair cutoff used by the potential
         int twojmax;     // larger values create more B components
         int switchflag;  // flag for smooth neighbour cutoff
         int bzeroflag;   // subtracts isolated atom bispectrum offsets
         int ncoeff;      // number of bispectrum coefficients (excluding the coefficient for j1=j2=j=0)

         void read_files();
         void read_coeff_file();
         void read_param_file();
         void build_material_map(const int num_materials);
         void print_info() const;
         int build_short_neighbour_list(const int atom,
                                        const std::vector<int>& neighbour_list_start_index,
                                        const std::vector<int>& neighbour_list_end_index,
                                        const std::vector<int>& type_array,
                                        const std::vector<int>& neighbour_list_array,
                                        const std::vector<double>& x_coord_array,
                                        const std::vector<double>& y_coord_array,
                                        const std::vector<double>& z_coord_array);
         void compute_bispectrum(const int start_index,
                                 const int end_index,
                                 const std::vector<int>& neighbour_list_start_index,
                                 const std::vector<int>& neighbour_list_end_index,
                                 const std::vector<int>& type_array,
                                 const std::vector<int>& neighbour_list_array,
                                 const std::vector<double>& x_coord_array,
                                 const std::vector<double>& y_coord_array,
                                 const std::vector<double>& z_coord_array,
                                 std::vector<double>& bispectrum);
         void compute_beta(const int start_index,
                           const int end_index,
                           const std::vector<int>& type_array,
                           std::vector<double>& beta) const;
         double compute_energy(const int atom,
                               const int local_atom,
                               const std::vector<int>& type_array,
                               const std::vector<double>& bispectrum) const;
         void print_debug_summary(const int call,
                                  const int start_index,
                                  const int end_index,
                                  const std::vector<int>& neighbour_counts,
                                  const std::vector<double>& initial_forces_x,
                                  const std::vector<double>& initial_forces_y,
                                  const std::vector<double>& initial_forces_z,
                                  const std::vector<double>& forces_array_x,
                                  const std::vector<double>& forces_array_y,
                                  const std::vector<double>& forces_array_z) const;
      };

   } // end of internal namespace

} // end of sld namespace

#endif // SLD_SNAP_H_
