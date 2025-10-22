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
// LS EDIT START
#include <limits>
#include <vector>
#include <fstream>
#include <iomanip>
// LS EDIT END

// Vampire headers
#include "atoms.hpp" // to be removed
#include "create.hpp" // to be removed

#include "anisotropy.hpp"
#include "errors.hpp"
#include "units.hpp"
#include "vio.hpp"

// anisotropy module headers
#include "internal.hpp"

namespace anisotropy{

   // LS EDIT START
   // Fn to output IDs of surface atoms detected 
   static inline void write_vampire_surface_csv(
      const std::vector<cs::catom_t>& catom_array)
   {
      static constexpr unsigned int FEA_ID = 0; // tetrahedral FeA material id
      static constexpr unsigned int FEB_ID = 1; // octahedral FeB material id

      const char* out_name = "vampire_surface.csv";
      std::ofstream ofs(out_name, std::ios::out | std::ios::trunc);
      if (!ofs) {
         zlog << zTs() << "WARNING: could not open " << out_name
               << " for writing." << std::endl;
         return;
      }

      // header
      ofs << "atom_id,type\n";

      // Only output rank local atoms 
      for (int a = 0; a < atoms::num_atoms; ++a) {
         if (catom_array[a].mpi_type == 2) continue; 

         if (!atoms::surface_array[a]) continue;     // only surface

         const unsigned int t = atoms::type_array[a];
         if (t == FEA_ID) {
               ofs << a << ",FeA\n";
         } else if (t == FEB_ID) {
               ofs << a << ",FeB\n";
         } else {
            
         }
      }

      ofs.close();
      zlog << zTs() << "Surface list written: " << out_name << std::endl;
   }
   
   // LS EDIT END

   //---------------------------------------------------------------------------
   // Function to identify less than fully coordinated atoms
   //---------------------------------------------------------------------------
   void identify_surface_atoms(std::vector<cs::catom_t> & catom_array, std::vector<std::vector <neighbours::neighbour_t> > & cneighbourlist){

      // initialise surface threshold if not overidden by input file
      if(internal::neel_anisotropy_threshold == 123456789) internal::neel_anisotropy_threshold = cs::unit_cell.surface_threshold;

      //-------------------------------------------------
      //	Optionally set up surface anisotropy
      //-------------------------------------------------

      // create temporary array for storing surface threshold
      std::vector<unsigned int> surface_anisotropy_threshold_array(atoms::num_atoms, internal::neel_anisotropy_threshold);
      // if using native (local) surface threshold then repopulate threshold array
      if(internal::native_neel_anisotropy_threshold){
         zlog << zTs() << "Identifying surface atoms using native (site dependent) threshold." << std::endl;
         for(int atom=0; atom < atoms::num_atoms; atom++){
            unsigned int atom_uc_id = catom_array.at(atom).uc_id;
            surface_anisotropy_threshold_array.at(atom) = cs::unit_cell.atom.at(atom_uc_id).ni;
         }
      }
      else zlog << zTs() << "Identifying surface atoms using global threshold value of " << internal::neel_anisotropy_threshold << std::endl;

      //--------------------------------------------------------------------------------------------
      // Determine nearest neighbour interactions from unit cell data for a single unit cell
      //--------------------------------------------------------------------------------------------

      // vector to store interactions within range
      std::vector<bool> nn_interaction(cs::unit_cell.bilinear.interaction.size(),false);

      // save nn_distance for performance
      const double rsq = internal::nearest_neighbour_distance * internal::nearest_neighbour_distance;

      // Get unit cell size
      const double ucdx = cs::unit_cell.dimensions[0];
      const double ucdy = cs::unit_cell.dimensions[1];
      const double ucdz = cs::unit_cell.dimensions[2];

      // loop over all interactions in unit cell
      for(unsigned int itr = 0; itr < cs::unit_cell.bilinear.interaction.size(); itr++){

         // get distance to neighbouring unit cell in unit cells
         double nndx = double(cs::unit_cell.bilinear.interaction[itr].dx);
         double nndy = double(cs::unit_cell.bilinear.interaction[itr].dy);
         double nndz = double(cs::unit_cell.bilinear.interaction[itr].dz);

         // load positions of i and j atoms to temporary coordinates and convert to angstroms
         double ix = (cs::unit_cell.atom[cs::unit_cell.bilinear.interaction[itr].i].x) * ucdx;
         double iy = (cs::unit_cell.atom[cs::unit_cell.bilinear.interaction[itr].i].y) * ucdy;
         double iz = (cs::unit_cell.atom[cs::unit_cell.bilinear.interaction[itr].i].z) * ucdz;
         double jx = (cs::unit_cell.atom[cs::unit_cell.bilinear.interaction[itr].j].x + nndx)*ucdx;
         double jy = (cs::unit_cell.atom[cs::unit_cell.bilinear.interaction[itr].j].y + nndy)*ucdy;
         double jz = (cs::unit_cell.atom[cs::unit_cell.bilinear.interaction[itr].j].z + nndz)*ucdz;

         // calculate reduced coordinates
         double dx = jx - ix;
         double dy = jy - iy;
         double dz = jz - iz;

         // calculate interaction range and check if less than nn distance
         const double range = (dx*dx + dy*dy + dz*dz);
         if(range <=rsq) nn_interaction[itr]=true;
      }

      //------------------------------------------------------------
      // Identify all nearest neighbour interactions in system
      //
      // Nearest neighbour list is a subset of full neighbour list,
      // and so everything is derived from that.
      //------------------------------------------------------------

      // vector to identify all nearest neighbour interactions
      std::vector <std::vector <bool> > nearest_neighbour_interactions_list(atoms::num_atoms);

      // loop over all atoms
      for(int atom=0; atom < atoms::num_atoms; atom++){

         // set all interactions for atom as non-nearest neighbour by default
         nearest_neighbour_interactions_list[atom].resize(cneighbourlist[atom].size(),false);

         // loop over all interactions for atom
         for(unsigned int nn=0;nn<cneighbourlist[atom].size();nn++){

            // get interaction type (same as unit cell interaction id)
            unsigned int id = cneighbourlist[atom][nn].i;

            // Ensure valid interaction id
            if(id>nn_interaction.size()){
               std::cout << "Error: invalid interaction id " << id << " is greater than number of interactions in unit cell " << nn_interaction.size() << ". Exiting" << std::endl;
               zlog << zTs() << "Error: invalid interaction id " << id << " is greater than number of interactions in unit cell " << nn_interaction.size() << ". Exiting" << std::endl;
               err::vexit();
            }

            // LS EDIT START keep but not using to determine periodicity
            nearest_neighbour_interactions_list[atom][nn]=nn_interaction.at(id);
            // LS EDIT END

         }
      }

      //----------------------------------------------------------------------------------------
      // Identify atoms with less than full nearest neighbour coordination
      //----------------------------------------------------------------------------------------

      // Track total number of surface atoms and total nearest neighbour interactions
      unsigned int num_surface_atoms=0;
      //unsigned int total_num_surface_nn=0;

      // Resize surface atoms mask and initialise to false
      atoms::surface_array.resize(atoms::num_atoms, false);

      // LS EDIT START
      static constexpr unsigned int FEA_ID = 0; // tetrahedral FeA material id
      static constexpr unsigned int FEB_ID = 1; // octahedral FeB material id
      static constexpr unsigned int O_ID   = 2; // oxygen material id
      static constexpr unsigned int THRESH_FEA = 4; // FeA coordination threshold (Fe–O only)
      static constexpr unsigned int THRESH_FEB = 6; // FeB coordination threshold (Fe–O only)
      // LS EDIT END

      // LS EDIT START Initial check of unit cell interactions for periodic hops and map to atoms - to identify which are from the xp xn yp yn surfaces
      const std::size_t num_uc_atoms = cs::unit_cell.atom.size();
      std::vector<unsigned int> uc_pbc_count(num_uc_atoms, 0);   // how many periodic interactions a site has
      std::vector<unsigned char> uc_has_pbc(num_uc_atoms, 0);    // flag if any

      // Bilinear
      for (std::size_t k = 0; k < cs::unit_cell.bilinear.interaction.size(); ++k) {
         const auto &I = cs::unit_cell.bilinear.interaction[k];
         if (I.dx != 0 || I.dy != 0 || I.dz != 0) {
            if (I.i >= 0 && (std::size_t)I.i < num_uc_atoms) { ++uc_pbc_count[I.i]; uc_has_pbc[I.i] = 1; }
            if (I.j >= 0 && (std::size_t)I.j < num_uc_atoms) { ++uc_pbc_count[I.j]; uc_has_pbc[I.j] = 1; }
         }
      }
      // Biquadratic (for completeness)
      for (std::size_t k = 0; k < cs::unit_cell.biquadratic.interaction.size(); ++k) {
         const auto &I = cs::unit_cell.biquadratic.interaction[k];
         if (I.dx != 0 || I.dy != 0 || I.dz != 0) {
            if (I.i >= 0 && (std::size_t)I.i < num_uc_atoms) { ++uc_pbc_count[I.i]; uc_has_pbc[I.i] = 1; }
            if (I.j >= 0 && (std::size_t)I.j < num_uc_atoms) { ++uc_pbc_count[I.j]; uc_has_pbc[I.j] = 1; }
         }
      }

      std::vector<unsigned char> atom_has_pbc(atoms::num_atoms, 0);
      std::vector<unsigned int>  atom_pbc_count(atoms::num_atoms, 0);
      unsigned int atoms_marked_pbc = 0;
      for (int a = 0; a < atoms::num_atoms; ++a) {
         const unsigned int ucid = catom_array[a].uc_id;
         if (ucid < uc_has_pbc.size()) {
            atom_has_pbc[a]  = uc_has_pbc[ucid];
            atom_pbc_count[a]= uc_pbc_count[ucid];
            if (atom_has_pbc[a]) ++atoms_marked_pbc;
         }
      }
      zlog << zTs() << "Pre-scan: " << atoms_marked_pbc
           << " atoms map to unitcellfile sites that have periodic hops (dx|dy|dz != 0)." << std::endl;
      // LS EDIT END

      // LS EDIT START store Fe–O coordination and debug fields
      std::vector<unsigned int> feo_coord(atoms::num_atoms, 0);
      std::vector<unsigned int> total_interactions_vec(atoms::num_atoms, 0);
      std::vector<unsigned char> excluded_due_to_pbc(atoms::num_atoms, 0);
      unsigned int undercoord_total = 0;
      unsigned int excluded_count   = 0;
      
      std::vector<unsigned int> pbc_any_nn(atoms::num_atoms, 0); 
      std::vector<unsigned int> pbc_feo_nn(atoms::num_atoms, 0); 
      // LS EDIT END

      // Single pass: classify surface iff (a) Fe–O undercoord and (b) not in precomputed PBC list
      for(int atom = 0; atom < atoms::num_atoms; atom++){

         // Check for local MPI atoms only
         if(catom_array[atom].mpi_type!=2){

            const unsigned int imat = atoms::type_array[atom];

            // Only classify Fe sites, O remain non-surface
            if (imat != FEA_ID && imat != FEB_ID){
               continue;
            }

            // Count Fe–O coordination using all entries in cneighbourlist since unitcell file already only includes nns
            unsigned int nFeO = 0;
            total_interactions_vec[atom] = static_cast<unsigned int>(cneighbourlist[atom].size());
            for (unsigned int nn = 0; nn < cneighbourlist[atom].size(); ++nn) {
               const unsigned int j_atom = cneighbourlist[atom][nn].nn;
               const unsigned int jmat   = atoms::type_array[j_atom];
               if (jmat == O_ID) ++nFeO;
            }
            feo_coord[atom] = nFeO;

            const unsigned int threshold = (imat == FEA_ID) ? THRESH_FEA : THRESH_FEB;
            const bool is_undercoord = (nFeO < threshold);

            if (is_undercoord) {
               ++undercoord_total;
               if (!atom_has_pbc[atom]) {
                  atoms::surface_array[atom] = true;
                  ++num_surface_atoms;
               } else {
                  atoms::surface_array[atom] = false;
                  excluded_due_to_pbc[atom] = 1;
                  ++excluded_count;
               }
            } else {
               atoms::surface_array[atom] = false;
            }

            // for debugging csv
            pbc_any_nn[atom] = atom_pbc_count[atom];
            pbc_feo_nn[atom] = 0; 
         } 
      }

      // Output statistics to log file
      zlog << zTs() << undercoord_total << " Fe atoms found under-coordinated (Fe–O)." << std::endl;
      zlog << zTs() << excluded_count   << " Fe atoms excluded due to periodic interactions." << std::endl;
      zlog << zTs() << num_surface_atoms<< " surface atoms after applying PBC exclusion." << std::endl;

      // LS EDIT START 
      {
         unsigned int nFeA = 0, nFeB = 0;
         for (int a = 0; a < atoms::num_atoms; ++a) if (atoms::surface_array[a]) {
            const unsigned int m = atoms::type_array[a];
            if (m == FEA_ID) ++nFeA; else if (m == FEB_ID) ++nFeB;
         }
         zlog << zTs() << num_surface_atoms << " surface Fe atoms found (FeA: "
              << nFeA << ", FeB: " << nFeB << ")" << std::endl;
      }
      // LS EDIT END

      // LS EDIT Optional csv showing surface atoms identified
      /*
      {
         const std::string csv_name = "surface_debug.csv";
         std::ofstream ofs(csv_name.c_str(), std::ios::out | std::ios::trunc);
         if (!ofs) {
            zlog << zTs() << "WARNING: could not open " << csv_name << " for writing." << std::endl;
         } else {
            // Columns: PBC count/flag, Fe–O coordination, total interactions, and exclusion reason
            ofs << "atom_id,uc_id,type,FeO_coord,threshold,"
                   "is_surface_after_pbc,n_pbc_uc,total_interactions,"
                   "excluded_due_to_pbc,x,y,z\n";
            ofs << std::setprecision(10);
            for (int a = 0; a < atoms::num_atoms; ++a) {
               const unsigned int t = atoms::type_array[a];
               if (t != FEA_ID && t != FEB_ID) continue; // only Fe sites
               const unsigned int thr = (t == FEA_ID) ? THRESH_FEA : THRESH_FEB;
               const bool is_surf_final = atoms::surface_array[a];
               const unsigned int ucid = catom_array[a].uc_id;
               const double x = catom_array[a].x;
               const double y = catom_array[a].y;
               const double z = catom_array[a].z;

               ofs << a << "," << ucid << "," << (t==FEA_ID?"FeA":"FeB") << ","
                   << feo_coord[a] << "," << thr << ","
                   << (is_surf_final?1:0) << ","
                   << atom_pbc_count[a] << ","
                   << total_interactions_vec[a] << ","
                   << static_cast<unsigned int>(excluded_due_to_pbc[a]) << ","
                   << x << "," << y << "," << z << "\n";
            }
            ofs.close();
            zlog << zTs() << "Surface debugging CSV written: " << csv_name << std::endl;
         }
      }
      */
      // LS EDIT END

      zlog << zTs() << "Surface atom identification complete." << std::endl;

      // LS EDIT START
      // Optionally write vampire_surface.csv
      {
         // Create a compact CSV with just the ids and Fe site type
         write_vampire_surface_csv(catom_array);
      }
      // LS EDIT END

      
      //----------------------------------------------------------------
      // If neel surface anisotropy is enabled, calculate necessary data
      //----------------------------------------------------------------
      if(internal::enable_neel_anisotropy){
         internal::initialise_neel_anisotropy_tensor(nearest_neighbour_interactions_list, cneighbourlist);
      }

      return;

   }

} // end of anisotropy namespace
