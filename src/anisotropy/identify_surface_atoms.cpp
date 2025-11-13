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

            // set mask to true or false for non-fully coordinated atoms in the bulk
            nearest_neighbour_interactions_list[atom][nn]=nn_interaction.at(id);

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


      // LS EDITS 
      static constexpr unsigned int FEA_ID = 0; // tetrahedral FeA material id
      static constexpr unsigned int FEB_ID = 1; // octahedral FeB material id
      static constexpr unsigned int O_ID = 2; // oxygen material id
      static constexpr unsigned int THRESH_FEA = 4; // FeA coordination threshold (Fe–O only)
      static constexpr unsigned int THRESH_FEB = 6; // FeB coordination threshold (Fe–O only)
      // END LS EDITS

      // Loop over all *local* atoms
      for(int atom = 0; atom < atoms::num_atoms; atom++){

         // Check for local MPI atoms only
         if(catom_array[atom].mpi_type!=2){

            // Initialise counter for number of nearest neighbour interactions
            unsigned int nnn_int=0;

            // Loop over all interactions to determine number of nearest neighbour interactions
            for(unsigned int nn = 0 ; nn < cneighbourlist[atom].size(); nn++){

               // If interaction is nn, increment counter
               if(nearest_neighbour_interactions_list[atom][nn]) nnn_int++;

            }

            // LS EDITS BELOW
            
            const unsigned int imat = atoms::type_array[atom];

            // Only classify Fe sites, Others (O) remain non-surface.
            if (imat != FEA_ID && imat != FEB_ID) continue;

            // If bulk Neel anisotropy is enabled, classify all Fe atoms as surface
            if (internal::enable_bulk_neel_anisotropy) {
               atoms::surface_array[atom] = true;
               ++num_surface_atoms;
               continue;
            }

            // Count only nearest neighbour oxygens
            unsigned int nnn_FeO = 0;
            for (unsigned int nn = 0; nn < cneighbourlist[atom].size(); ++nn) {
               // only consider interactions labelled as nearest neighbour
               if (!nearest_neighbour_interactions_list[atom][nn]) continue;

               const unsigned int j_atom = cneighbourlist[atom][nn].nn; // neighbour atom index
               const unsigned int jmat = atoms::type_array[j_atom]; // neighbour material id
               if (jmat == O_ID) ++nnn_FeO; // count Fe–O only
            }

            // Choose Fe site specific threshold and override any global threshold
            unsigned int threshold = (imat == FEA_ID) ? THRESH_FEA : THRESH_FEB;

            if (nnn_FeO < threshold) {
               atoms::surface_array[atom] = true;
               ++num_surface_atoms;
            }
         }
         // END LS EDITS


         
      }

      // Output statistics to log file
      // LS EDITS(counts by FeA/FeB)
      {
         unsigned int nFeA = 0, nFeB = 0;
         for (int a = 0; a < atoms::num_atoms; ++a) {
            if (atoms::surface_array[a]) {
               const unsigned int m = atoms::type_array[a];
               if (m == FEA_ID) ++nFeA; 
               else if (m == FEB_ID) ++nFeB;
            }
         }
         if (internal::enable_bulk_neel_anisotropy) {
            zlog << zTs() << "Bulk Neel anisotropy enabled: All " << num_surface_atoms 
                 << " Fe atoms classified as surface (FeA: " << nFeA << ", FeB: " << nFeB << ")" << std::endl;
         }
         else {
            zlog << zTs() << num_surface_atoms << " surface Fe atoms found (FeA: "
                 << nFeA << ", FeB: " << nFeB << ")" << std::endl;
         }
      }
      // END LS EDITS
      
      //----------------------------------------------------------------
      // If neel surface anisotropy is enabled, calculate necessary data
      //----------------------------------------------------------------
      if(internal::enable_neel_anisotropy){
         internal::initialise_neel_anisotropy_tensor(nearest_neighbour_interactions_list, cneighbourlist);
      }

      return;

   }

} // end of anisotropy namespace
