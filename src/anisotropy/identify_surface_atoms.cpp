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
#include <limits>
#include <vector>
#include <fstream>
#include <iomanip>

// Vampire headers
#include "atoms.hpp"   // to be removed
#include "create.hpp"  // to be removed

#include "anisotropy.hpp"
#include "errors.hpp"
#include "units.hpp"
#include "vio.hpp"

// anisotropy module headers
#include "internal.hpp"

namespace anisotropy{

   // ---------------- constants ----------------
   static constexpr unsigned int FEA_ID = 0; // FeA (tet)
   static constexpr unsigned int FEB_ID = 1; // FeB (oct)
   static constexpr unsigned int O_ID   = 2; // O

   // Surface thresholds by site:
   static constexpr unsigned int THRESH_FEA = 4; // FeA (tet)
   static constexpr unsigned int THRESH_FEB = 6; // FeB (oct)

   // ---------------- helpers ----------------
   // compact surface atom list
   static inline void write_vampire_surface_csv(
      const std::vector<cs::catom_t>& catom_array)
   {
      const char* out_name = "vampire_surface.csv";
      std::ofstream ofs(out_name, std::ios::out | std::ios::trunc);
      if (!ofs) {
         zlog << zTs() << "WARNING: could not open " << out_name
              << " for writing." << std::endl;
         return;
      }

      ofs << "atom_id,type\n";
      for (int a = 0; a < atoms::num_atoms; ++a) {
         if (catom_array[a].mpi_type == 2) continue;  
         if (!atoms::surface_array[a])    continue;   // only surface

         const unsigned int t = atoms::type_array[a];
         if (t == FEA_ID)      ofs << a << ",FeA\n";
         else if (t == FEB_ID) ofs << a << ",FeB\n";
      }

      ofs.close();
      zlog << zTs() << "Surface list written: " << out_name << std::endl;
   }

//---------------------------------------------------------------------------
// Identify surface atoms of a slab of magnetite 
//  FeB is surface if Fe–O interactions < 6
//  FeA is surface if Fe–O interactions < 4
//  any Fe that has any periodic interaction is not allowed to be surface.
//---------------------------------------------------------------------------
   void identify_surface_atoms(std::vector<cs::catom_t> & catom_array,
                            std::vector<std::vector <neighbours::neighbour_t> > & cneighbourlist)
   {
      // initialise surface threshold if not overridden by input file
      if(internal::neel_anisotropy_threshold == 123456789)
         internal::neel_anisotropy_threshold = cs::unit_cell.surface_threshold;

      // create temporary array for storing surface threshold
      std::vector<unsigned int> surface_anisotropy_threshold_array(
         atoms::num_atoms, internal::neel_anisotropy_threshold);

      // if using native threshold then repopulate threshold array
      if(internal::native_neel_anisotropy_threshold){
         zlog << zTs() << "Identifying surface atoms using native (site dependent) threshold." << std::endl;
         for(int atom=0; atom < atoms::num_atoms; atom++){
            unsigned int atom_uc_id = catom_array.at(atom).uc_id;
            surface_anisotropy_threshold_array.at(atom) = cs::unit_cell.atom.at(atom_uc_id).ni;
         }
      }
      else {
         zlog << zTs() << "Identifying surface atoms using global threshold value of "
            << internal::neel_anisotropy_threshold << std::endl;
      }

      //--------------------------------------------------------------------------------------------
      // Build nearest_neighbour_interactions_list shape (kept for downstream tensor init)
      //--------------------------------------------------------------------------------------------
      std::vector <std::vector <bool> > nearest_neighbour_interactions_list(atoms::num_atoms);
      for(int atom=0; atom < atoms::num_atoms; atom++){
         nearest_neighbour_interactions_list[atom].resize(cneighbourlist[atom].size(), false);
         for(unsigned int nn=0; nn<cneighbourlist[atom].size(); nn++){
            nearest_neighbour_interactions_list[atom][nn] = true; // NN by construction in your input
         }
      }

      //----------------------------------------------------------------------------------------
      // 1) Build periodic mask ... mark any Fe site with any periodic hop
      //----------------------------------------------------------------------------------------
      const std::size_t num_uc_atoms = cs::unit_cell.atom.size();
      std::vector<unsigned int>  uc_pbc_count(num_uc_atoms, 0);
      std::vector<unsigned char> uc_has_pbc (num_uc_atoms, 0);

      auto mark_fe_if_periodic = [&](unsigned int uidx, bool isFe, int dx, int dy, int dz){
         if(!isFe) return;
         if (dx != 0 || dy != 0 || dz != 0) {
            ++uc_pbc_count[uidx];
            uc_has_pbc[uidx] = 1;
         }
      };

      // Bilinear
      for (std::size_t k = 0; k < cs::unit_cell.bilinear.interaction.size(); ++k) {
         const auto &I  = cs::unit_cell.bilinear.interaction[k];
         const unsigned int ui = static_cast<unsigned int>(I.i);
         const unsigned int uj = static_cast<unsigned int>(I.j);
         if (ui >= num_uc_atoms || uj >= num_uc_atoms) continue;

         const unsigned int ti = cs::unit_cell.atom[ui].mat;
         const unsigned int tj = cs::unit_cell.atom[uj].mat;
         const bool uiIsFe = (ti == FEA_ID || ti == FEB_ID);
         const bool ujIsFe = (tj == FEB_ID || tj == FEA_ID);

         mark_fe_if_periodic(ui, uiIsFe, I.dx, I.dy, I.dz);
         mark_fe_if_periodic(uj, ujIsFe, I.dx, I.dy, I.dz);
      }

      // Biquadratic
      for (std::size_t k = 0; k < cs::unit_cell.biquadratic.interaction.size(); ++k) {
         const auto &I  = cs::unit_cell.biquadratic.interaction[k];
         const unsigned int ui = static_cast<unsigned int>(I.i);
         const unsigned int uj = static_cast<unsigned int>(I.j);
         if (ui >= num_uc_atoms || uj >= num_uc_atoms) continue;

         const unsigned int ti = cs::unit_cell.atom[ui].mat;
         const unsigned int tj = cs::unit_cell.atom[uj].mat;
         const bool uiIsFe = (ti == FEA_ID || ti == FEB_ID);
         const bool ujIsFe = (tj == FEB_ID || tj == FEA_ID);

         mark_fe_if_periodic(ui, uiIsFe, I.dx, I.dy, I.dz);
         mark_fe_if_periodic(uj, ujIsFe, I.dx, I.dy, I.dz);
      }

      std::vector<unsigned char> atom_has_pbc(atoms::num_atoms, 0);
      std::vector<unsigned int>  atom_pbc_count(atoms::num_atoms, 0);
      unsigned int atoms_marked_pbc = 0;
      for (int a = 0; a < atoms::num_atoms; ++a) {
         const unsigned int ucid = catom_array[a].uc_id;
         if (ucid < uc_has_pbc.size()) {
            atom_has_pbc[a]   = uc_has_pbc[ucid];
            atom_pbc_count[a] = uc_pbc_count[ucid];
            if (atom_has_pbc[a]) ++atoms_marked_pbc;
         }
      }
      zlog << zTs() << "Pre-scan: " << atoms_marked_pbc
         << " Fe atoms map to UC sites with at least one periodic interaction (dx|dy|dz != 0)."
         << std::endl;

      //----------------------------------------------------------------------------------------
      // 2) Fe–O coordination per site from unit_cell interactions 
      //----------------------------------------------------------------------------------------
      std::vector<unsigned int> feo_coord_uc_raw(num_uc_atoms, 0);

      for (unsigned int k = 0; k < cs::unit_cell.bilinear.interaction.size(); ++k) {
         const auto &I  = cs::unit_cell.bilinear.interaction[k];
         const unsigned int ui = static_cast<unsigned int>(I.i);
         const unsigned int uj = static_cast<unsigned int>(I.j);
         if (ui >= num_uc_atoms || uj >= num_uc_atoms) continue;

         const unsigned int ti = cs::unit_cell.atom[ui].mat;
         const unsigned int tj = cs::unit_cell.atom[uj].mat;

         const bool iFe = (ti == FEA_ID || ti == FEB_ID);
         const bool jFe = (tj == FEA_ID || tj == FEB_ID);
         const bool iO  = (ti == O_ID);
         const bool jO  = (tj == O_ID);

         if (iFe && jO) ++feo_coord_uc_raw[ui]; // Fe(i) -> O(j)
         if (jFe && iO) ++feo_coord_uc_raw[uj]; // O(i) -> Fe(j)
      }

      //----------------------------------------------------------------------------------------
      // 3) Classify surface atoms using halved Fe–O counts + PBC exclusion
      //----------------------------------------------------------------------------------------
      atoms::surface_array.resize(atoms::num_atoms, false);

      std::vector<unsigned int>  feo_coord_raw(atoms::num_atoms, 0);     // debug ... raw double counted
      std::vector<unsigned int>  feo_coord_half(atoms::num_atoms, 0);    // used for decision
      std::vector<unsigned int>  total_interactions_vec(atoms::num_atoms, 0); // debug
      std::vector<unsigned char> excluded_due_to_pbc(atoms::num_atoms, 0);

      unsigned int undercoord_total = 0;
      unsigned int excluded_count   = 0;
      unsigned int num_surface_atoms= 0;

      for(int atom = 0; atom < atoms::num_atoms; atom++){
         if(catom_array[atom].mpi_type == 2) continue;

         const unsigned int imat = atoms::type_array[atom];
         if (imat != FEA_ID && imat != FEB_ID) continue; // Fe only

         const unsigned int ucid = catom_array[atom].uc_id;

         const unsigned int nFeO_raw  = (ucid < feo_coord_uc_raw.size()) ? feo_coord_uc_raw[ucid] : 0u;
         const unsigned int nFeO_half = nFeO_raw / 2u; // remove i-j / j-i duplication

         feo_coord_raw[atom]  = nFeO_raw;
         feo_coord_half[atom] = nFeO_half;

         // debug how many neighbours in cneighbourlist
         total_interactions_vec[atom] = static_cast<unsigned int>(cneighbourlist[atom].size());

         const unsigned int threshold = (imat == FEA_ID) ? THRESH_FEA : THRESH_FEB;
         const bool is_undercoord = (nFeO_half < threshold);

         if (is_undercoord) {
            ++undercoord_total;
            // PBC exclusion if this Fe has periodic interactions, it's not allowed to be surface
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
      }

      // Stats
      zlog << zTs() << undercoord_total << " Fe atoms found under-coordinated by Fe–O counts." << std::endl;
      zlog << zTs() << excluded_count   << " under-coordinated Fe excluded due to periodic interactions." << std::endl;
      zlog << zTs() << num_surface_atoms<< " surface Fe atoms after applying PBC exclusion." << std::endl;

      {
         unsigned int nFeA = 0, nFeB = 0;
         for (int a = 0; a < atoms::num_atoms; ++a) if (atoms::surface_array[a]) {
            const unsigned int m = atoms::type_array[a];
            if (m == FEA_ID) ++nFeA; else if (m == FEB_ID) ++nFeB;
         }
         zlog << zTs() << num_surface_atoms << " surface Fe atoms found (FeA: "
            << nFeA << ", FeB: " << nFeB << ")" << std::endl;
      }

      // -------- Debugging CSVs --------
      // {
      //    const std::string csv_name = "surface_debug.csv";
      //    std::ofstream ofs(csv_name.c_str(), std::ios::out | std::ios::trunc);
      //    if (!ofs) {
      //       zlog << zTs() << "WARNING: could not open " << csv_name << " for writing." << std::endl;
      //    } else {
      //       ofs << "atom_id,uc_id,type,FeO_coord_raw,FeO_coord_used,threshold,"
      //             "is_surface_after_pbc,n_pbc_uc,total_interactions,"
      //             "excluded_due_to_pbc,x,y,z\n";
      //       ofs << std::setprecision(10);
      //       for (int a = 0; a < atoms::num_atoms; ++a) {
      //          const unsigned int t = atoms::type_array[a];
      //          if (t != FEA_ID && t != FEB_ID) continue; // Fe only
      //          const unsigned int thr = (t == FEA_ID) ? THRESH_FEA : THRESH_FEB;
      //          const bool is_surf_final = atoms::surface_array[a];
      //          const unsigned int ucid = catom_array[a].uc_id;
      //          const double x = catom_array[a].x;
      //          const double y = catom_array[a].y;
      //          const double z = catom_array[a].z;

      //          ofs << a << "," << ucid << "," << (t==FEA_ID?"FeA":"FeB") << ","
      //             << feo_coord_raw[a]  << "," << feo_coord_half[a] << ","
      //             << thr << ","
      //             << (is_surf_final?1:0) << ","
      //             << atom_pbc_count[a] << ","
      //             << total_interactions_vec[a] << ","
      //             << static_cast<unsigned int>(excluded_due_to_pbc[a]) << ","
      //             << x << "," << y << "," << z << "\n";
      //       }
      //       ofs.close();
      //       zlog << zTs() << "Surface debugging CSV written: " << csv_name << std::endl;
      //    }
      // }

      // compact surface list (FeA/FeB only)
      write_vampire_surface_csv(catom_array);

      // If neel surface anisotropy is enabled, calculate necessary data
      if(internal::enable_neel_anisotropy){
         internal::initialise_neel_anisotropy_tensor(
            nearest_neighbour_interactions_list, cneighbourlist);
      }

      return;
   }

} // end of anisotropy namespace
