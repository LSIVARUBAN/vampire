//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Sam Westmoreland and Richard Evans 2017. All rights reserved.
//   (c) Lennon Sivaruban 2026. All rights reserved.
//
//   Email: richard.evans@york.ac.uk
//          lennon.sivaruban@manchester.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <string>
#include <sstream>

// Vampire headers
#include "atoms.hpp"
#include "create.hpp"
#include "errors.hpp"
#include "material.hpp"
#include "vio.hpp"

// create module headers
#include "internal.hpp"

namespace create{
namespace internal{

   //---------------------------------------------------------------------------
   // Function to identify less than fully coordinated (surface) atoms in the
   // initial structure, saved in atoms::surface_array, and the nearest neighbour
   // interactions of each atom (used for Neel anisotropy)
   //---------------------------------------------------------------------------
   void identify_surface_atoms(std::vector<cs::catom_t> & catom_array,
                               std::vector<std::vector <neighbours::neighbour_t> > & cneighbourlist,
                               std::vector <std::vector <bool> >& nearest_neighbour_interactions_list){

      // initialise surface threshold if not overidden by input file
      if(internal::surface_threshold == 123456789) internal::surface_threshold = cs::unit_cell.surface_threshold;

      // use surface flags from the unit cell file if present
      const bool use_unit_cell_surface_flags = cs::unit_cell.surface_flags_present;

      // create temporary array for storing surface threshold
      std::vector<unsigned int> surface_threshold_array(atoms::num_atoms, internal::surface_threshold);
      if(use_unit_cell_surface_flags) zlog << zTs() << "Identifying surface atoms using unit cell surface flags." << std::endl;
      else zlog << zTs() << "Identifying surface atoms using global threshold value of " << internal::surface_threshold << std::endl;

      //--------------------------------------------------------------------------------------------
      // Determine nearest neighbour interactions from unit cell data for a single unit cell
      //--------------------------------------------------------------------------------------------

      // vector to store interactions within range
      std::vector<bool> nn_interaction(cs::unit_cell.bilinear.interaction.size(),false);

      // nearest neighbour range (A), by default the nearest neighbour radius of a built-in unit cell
      // (so that longer range exchange interactions are not counted) or all interactions in a unit cell file
      double nn_range = internal::surface_identification_range;
      if(nn_range < 0.0) nn_range = (cs::unit_cell.cutoff_radius > 0.0) ? 1.000001*cs::unit_cell.cutoff_radius*cs::unit_cell.dimensions[0] : 1.0e9;
      zlog << zTs() << "Counting nearest neighbours within " << nn_range << " A for surface atom identification." << std::endl;

      // save nn_distance for performance
      const double rsq = nn_range * nn_range;

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
      nearest_neighbour_interactions_list.resize(atoms::num_atoms);

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

      // Loop over all *local* atoms
      for(int atom = 0; atom < atoms::num_atoms; atom++){

         // Check for local MPI atoms only
         if(catom_array[atom].mpi_type!=2){

            // surface atoms defined by the unit cell file flags
            if(use_unit_cell_surface_flags){
               if(cs::unit_cell.atom[catom_array[atom].uc_id].is_surface){
                  atoms::surface_array[atom]=true;
                  num_surface_atoms++;
               }
               continue;
            }

            // Initialise counter for number of nearest neighbour interactions
            unsigned int nnn_int=0;

            // Loop over all interactions to determine number of nearest neighbour interactions
            for(unsigned int nn = 0 ; nn < cneighbourlist[atom].size(); nn++){

               // If interaction is nn, increment counter
               if(nearest_neighbour_interactions_list[atom][nn]) nnn_int++;

            }

            // check for atoms with < threshold number of nearest neighbours
            if(nnn_int<surface_threshold_array.at(atom)){
               atoms::surface_array[atom]=true;
               num_surface_atoms++;
               //total_num_surface_nn+=nnn_int;
            }
         }
      }

      // Output statistics to log file
      zlog << zTs() << num_surface_atoms << " surface atoms found." << std::endl;

      // Output number of surface atoms for each material to log file
      std::vector<unsigned int> num_surface_atoms_per_material(mp::num_materials, 0);
      for(int atom = 0; atom < atoms::num_atoms; atom++){
         if(atoms::surface_array[atom]) num_surface_atoms_per_material[atoms::type_array[atom]]++;
      }
      zlog << zTs() << "Number of surface atoms for each material:";
      for(int mat = 0; mat < mp::num_materials; mat++) zlog << " " << num_surface_atoms_per_material[mat];
      zlog << std::endl;

      return;

   }

} // end of internal namespace
} // end of create namespace
