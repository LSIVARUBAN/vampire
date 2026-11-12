//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Lennon Sivaruban 2026. All rights reserved.
//
//   Email: lennon.sivaruban@manchester.ac.uk
//
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <iomanip>
#include <sstream>

// Vampire headers
#include "atoms.hpp"
#include "config.hpp"
#include "material.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "vmpi.hpp"

// config module headers
#include "internal.hpp"

namespace config{

namespace internal{

//------------------------------------------------------------------------------
// Function to set up output of the surface atoms (atoms::surface_array) among
// the atoms to be outputted. Surface data are gathered on the root process and
// written to a single file in the same format (and with the same meta files) as
// the corresponding atom files.
//------------------------------------------------------------------------------
void initialize_surface(){

   // list of local surface atoms to output
   for(uint64_t id = 0; id < local_output_atom_list.size(); id++){
      if(atoms::surface_array[local_output_atom_list[id]]) local_surface_atom_list.push_back(local_output_atom_list[id]);
   }

   // number of surface atoms on each processor and their offsets in the gathered data on the root process
   int local_atoms = local_surface_atom_list.size();
   surface_recv_counts.assign(vmpi::num_processors, 0);
   surface_displacements.assign(vmpi::num_processors, 0);
   #ifdef MPICF
      MPI_Gather(&local_atoms, 1, MPI_INT, &surface_recv_counts[0], 1, MPI_INT, 0, MPI_COMM_WORLD);
   #else
      surface_recv_counts[0] = local_atoms;
   #endif
   total_surface_atoms = 0;
   for(int p = 0; p < vmpi::num_processors; p++){
      surface_displacements[p] = total_surface_atoms;
      total_surface_atoms += surface_recv_counts[p];
   }

   zlog << zTs() << "Outputting configurations of " << total_surface_atoms << " surface atoms." << std::endl;

   return;

}

//------------------------------------------------------------------------------
// Function to gather data of the local surface atoms on the root process, with
// num_values values per atom (3 for vectors)
//------------------------------------------------------------------------------
template <typename T> void gather_surface_data(const std::vector<T>& local, std::vector<T>& gathered, const int num_values){

   #ifdef MPICF
      std::vector<int> counts(surface_recv_counts.size());
      std::vector<int> displacements(surface_displacements.size());
      for(unsigned int p = 0; p < counts.size(); p++){
         counts[p] = num_values * surface_recv_counts[p];
         displacements[p] = num_values * surface_displacements[p];
      }
      if(vmpi::my_rank == 0) gathered.resize(num_values * total_surface_atoms);
      const MPI_Datatype type = (sizeof(T) == sizeof(double)) ? MPI_DOUBLE : MPI_INT;
      MPI_Gatherv(local.data(), local.size(), type, gathered.data(), counts.data(), displacements.data(), type, 0, MPI_COMM_WORLD);
   #else
      gathered = local;
   #endif

   return;

}

//------------------------------------------------------------------------------
// Function to gather the 3-vector data (x,y,z) of the surface atoms on the root process
//------------------------------------------------------------------------------
void gather_surface_vectors(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z, std::vector<double>& gathered){

   std::vector<double> local(3 * local_surface_atom_list.size());
   copy_data_to_buffer(x, y, z, local_surface_atom_list, local);
   gather_surface_data(local, gathered, 3);

   return;

}

//------------------------------------------------------------------------------
// Function to write the meta file of a configuration of the surface atoms
// (spins-surface-XXXXXXXX.meta, same format as spins-XXXXXXXX.meta)
//------------------------------------------------------------------------------
void write_surface_meta(){

   const std::vector<double> magnetisation = stats::system_magnetization.get_magnetization();
   const double real_time = double(sim::time) * mp::dt_SI;
   write_meta(real_time, sim::temperature, sim::H_vec[0], sim::H_vec[1], sim::H_vec[2], sim::H_applied, magnetisation[0], magnetisation[1], magnetisation[2], true);

   return;

}

//------------------------------------------------------------------------------
// Function to output the material, category and coordinates of the surface atoms
// (atoms-coords-surface.data and .meta, same format as atoms-coords.data and .meta)
//------------------------------------------------------------------------------
void atoms_coords_surface(){

   // material and category of local surface atoms
   std::vector<int> local_type(local_surface_atom_list.size());
   std::vector<int> local_category(local_surface_atom_list.size());
   for(unsigned int id = 0; id < local_surface_atom_list.size(); id++){
      local_type[id] = atoms::type_array[local_surface_atom_list[id]];
      local_category[id] = atoms::category_array[local_surface_atom_list[id]];
   }

   std::vector<double> coords;
   std::vector<int> type;
   std::vector<int> category;
   gather_surface_vectors(atoms::x_coord_array, atoms::y_coord_array, atoms::z_coord_array, coords);
   gather_surface_data(local_type, type, 1);
   gather_surface_data(local_category, category, 1);

   if(vmpi::my_rank == 0){
      write_coordinate_meta(true);
      write_coord_data("atoms-coords-surface.data", coords, type, category);
   }

   return;

}

//------------------------------------------------------------------------------
// Function to output the spins of the surface atoms
// (spins-surface-XXXXXXXX.data and .meta, same format as spins-XXXXXXXX.data and .meta)
//------------------------------------------------------------------------------
void atoms_surface(){

   std::vector<double> spins;
   gather_surface_vectors(atoms::x_spin_array, atoms::y_spin_array, atoms::z_spin_array, spins);

   std::stringstream filename;
   filename << "spins-surface-" << std::setfill('0') << std::setw(8) << sim::output_atoms_file_counter << ".data";
   if(vmpi::my_rank == 0){
      write_surface_meta();
      write_data(filename.str(), spins);
   }

   return;

}

//------------------------------------------------------------------------------
// Function to output the dynamic (SLD) positions of the surface atoms
// (atoms-coords-sld-surface-XXXXXXXX.data, same format as atoms-coords-sld-XXXXXXXX.data,
// with the configuration meta file as for atoms-coords-sld)
//------------------------------------------------------------------------------
void atoms_coords_sld_surface(){

   std::vector<double> coords;
   gather_surface_vectors(atoms::x_coord_array, atoms::y_coord_array, atoms::z_coord_array, coords);

   std::stringstream filename;
   filename << "atoms-coords-sld-surface-" << std::setfill('0') << std::setw(8) << sim::output_atoms_file_counter << ".data";
   if(vmpi::my_rank == 0){
      write_surface_meta();
      write_data(filename.str(), coords);
   }

   return;

}

} // end of internal namespace

} // end of config namespace
