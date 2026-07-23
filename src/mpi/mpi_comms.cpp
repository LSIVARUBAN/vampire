//-----------------------------------------------------------------------------
//
//  Vampire - A code for atomistic simulation of magnetic materials
//
//  Copyright (C) 2009-2012 R.F.L.Evans
//
//  Email:richard.evans@york.ac.uk
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful, but
//  WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
//  General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software Foundation,
//  Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
//
// ----------------------------------------------------------------------------
//
//====================================================================================
//
//														mpi_comms
//
//						Functions for performing halo swap during parallel execution
//
//										Version 0.1 R Evans 14/09/2009
//
//====================================================================================
//
//		Locally allocated variables:
//
//=====================================================================================

#include "atoms.hpp"
#include "errors.hpp"
#include "vmpi.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>

namespace vmpi{

void mpi_init_halo_swap(){
	//====================================================================================
	//
	///												mpi_init_halo_swap
	///
	///									Initiates halo swap for spin data
	///
	///										Version 1.0 R Evans 16/09/2009
	//
	//====================================================================================
	//
	//		Locally allocated variables:
	//
	//====================================================================================
	//using namespace mpi_comms;
#ifdef MPICF
	//----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
	if(err::check==true){
		std::cout << "mpi_init_halo_swap has been called" << "\t";
		std::cout << vmpi::my_rank << std::endl;
	}

	//----------------------------------------------------------
	// Pack spins for sending
	//----------------------------------------------------------

	//std::cout << "rank: " << vmpi::my_rank << " num_halo_swaps: " << vmpi::send_atom_translation_array.size();
	//std::cout << "\t" << "num_boundary_swaps: " << vmpi::recv_atom_translation_array.size() << std::endl;

	for(unsigned int i=0;i<vmpi::send_atom_translation_array.size();i++){
		int atom = vmpi::send_atom_translation_array[i];
		vmpi::send_spin_data_array[3*i+0] = atoms::x_spin_array[atom];
		vmpi::send_spin_data_array[3*i+1] = atoms::y_spin_array[atom];
		vmpi::send_spin_data_array[3*i+2] = atoms::z_spin_array[atom];
		//std::cout << vmpi::my_rank << "\t" << i << "\t";
		//std::cout << atoms::x_spin_array[vmpi::send_atom_translation_array[i]] << "\t";
		//std::cout << atoms::y_spin_array[vmpi::send_atom_translation_array[i]] << "\t";
		//std::cout << atoms::z_spin_array[vmpi::send_atom_translation_array[i]] << std::endl;
	}

	//----------------------------------------------------------
	// Send spin data array and post receives
	//----------------------------------------------------------

	vmpi::requests.resize(0);
	MPI_Request req = 0;

	for (int p=0;p<vmpi::num_processors;p++){
		if(vmpi::send_num_array[p]!=0){
			int num_pts = 3*vmpi::send_num_array[p];
			int si = 3*vmpi::send_start_index_array[p];
			vmpi::requests.push_back(req);
			MPI_Isend(&vmpi::send_spin_data_array[si],num_pts,MPI_DOUBLE,p,48, MPI_COMM_WORLD, &vmpi::requests.back());
		}
		if(vmpi::recv_num_array[p]!=0){
			int num_pts = 3*vmpi::recv_num_array[p];
			int si = 3*vmpi::recv_start_index_array[p];
			vmpi::requests.push_back(req);
			MPI_Irecv(&vmpi::recv_spin_data_array[si],num_pts,MPI_DOUBLE,p,48, MPI_COMM_WORLD, &vmpi::requests.back());
		}
	}

   #endif
	//----------------------------------------------------------
	// Return
	//----------------------------------------------------------

	return;

}

void mpi_complete_halo_swap(){
	//====================================================================================
	///
	///												mpi_complete_halo_swap
	///
	///									Completes halo swap for spin data
	///
	///										Version 1.0 R Evans 16/09/2009
	//
	//====================================================================================

   #ifdef MPICF

	// check calling of routine if error checking is activated
	if(err::check==true){
		std::cout << "mpi_complete_halo_swap has been called" << "\t";
		std::cout << vmpi::my_rank << std::endl;
	}

	// Swap timers compute -> wait
	vmpi::TotalComputeTime+=vmpi::SwapTimer(vmpi::ComputeTime, vmpi::WaitTime);

	// Wait for all comms to complete
	vmpi::stati.resize(vmpi::requests.size());
	MPI_Waitall(vmpi::requests.size(),&vmpi::requests[0],&vmpi::stati[0]);

	// Swap timers wait -> compute
	vmpi::TotalWaitTime+=vmpi::SwapTimer(vmpi::WaitTime, vmpi::ComputeTime);

	// Unpack received spins
	for(unsigned int i=0;i<vmpi::recv_atom_translation_array.size();i++){
		atoms::x_spin_array[vmpi::recv_atom_translation_array[i]] = vmpi::recv_spin_data_array[3*i+0];
		atoms::y_spin_array[vmpi::recv_atom_translation_array[i]] = vmpi::recv_spin_data_array[3*i+1];
		atoms::z_spin_array[vmpi::recv_atom_translation_array[i]] = vmpi::recv_spin_data_array[3*i+2];
		//std::cout << vmpi::my_rank << "\t" << i << "\t" << vmpi::recv_atom_translation_array[i] << "\t";
		//std::cout << atoms::x_spin_array[vmpi::recv_atom_translation_array[i]] << "\t";
		//std::cout << atoms::y_spin_array[vmpi::recv_atom_translation_array[i]] << "\t";
		//std::cout << atoms::z_spin_array[vmpi::recv_atom_translation_array[i]] << std::endl;
	}

   #endif

	return;

}

void mpi_init_halo_swap_coords(){
#ifdef MPICF
	//----------------------------------------------------------
	// check calling of routine if error checking is activated
	//----------------------------------------------------------
	if(err::check==true){
		std::cout << "mpi_init_halo_swap has been called" << "\t";
		std::cout << vmpi::my_rank << std::endl;
	}



	//----------------------------------------------------------
	// Pack spins for sending
	//----------------------------------------------------------

	//std::cout << "rank: " << vmpi::my_rank << " num_halo_swaps: " << vmpi::send_atom_translation_array.size();
	//std::cout << "\t" << "num_boundary_swaps: " << vmpi::recv_atom_translation_array.size() << std::endl;

	for(unsigned int i=0;i<vmpi::send_atom_translation_array.size();i++){
		int atom = vmpi::send_atom_translation_array[i];
		vmpi::send_coord_data_array[3*i+0] = atoms::x_coord_array[atom];
		vmpi::send_coord_data_array[3*i+1] = atoms::y_coord_array[atom];
		vmpi::send_coord_data_array[3*i+2] = atoms::z_coord_array[atom];


	}

	//----------------------------------------------------------
	// Send spin data array and post receives
	//----------------------------------------------------------

	vmpi::requests.resize(0);
	MPI_Request req = 0;

	for (int p=0;p<vmpi::num_processors;p++){
		if(vmpi::send_num_array[p]!=0){
			int num_pts = 3*vmpi::send_num_array[p];
			int si = 3*vmpi::send_start_index_array[p];
			vmpi::requests.push_back(req);
			MPI_Isend(&vmpi::send_coord_data_array[si],num_pts,MPI_DOUBLE,p,48, MPI_COMM_WORLD, &vmpi::requests.back());
		}
		if(vmpi::recv_num_array[p]!=0){
			int num_pts = 3*vmpi::recv_num_array[p];
			int si = 3*vmpi::recv_start_index_array[p];
			vmpi::requests.push_back(req);
			MPI_Irecv(&vmpi::recv_coord_data_array[si],num_pts,MPI_DOUBLE,p,48, MPI_COMM_WORLD, &vmpi::requests.back());
		}
	}

   #endif
	//----------------------------------------------------------
	// Return
	//----------------------------------------------------------

	return;

}

void mpi_complete_halo_swap_coords(){
   #ifdef MPICF

	// check calling of routine if error checking is activated
	if(err::check==true){
		std::cout << "mpi_complete_halo_swap has been called" << "\t";
		std::cout << vmpi::my_rank << std::endl;
	}


	// Swap timers compute -> wait
	vmpi::TotalComputeTime+=vmpi::SwapTimer(vmpi::ComputeTime, vmpi::WaitTime);

	// Wait for all comms to complete
	vmpi::stati.resize(vmpi::requests.size());
	MPI_Waitall(vmpi::requests.size(),&vmpi::requests[0],&vmpi::stati[0]);

	// Swap timers wait -> compute
	vmpi::TotalWaitTime+=vmpi::SwapTimer(vmpi::WaitTime, vmpi::ComputeTime);

	// Unpack received spins
	for(unsigned int i=0;i<vmpi::recv_atom_translation_array.size();i++){
		atoms::x_coord_array[vmpi::recv_atom_translation_array[i]] = vmpi::recv_coord_data_array[3*i+0];
		atoms::y_coord_array[vmpi::recv_atom_translation_array[i]] = vmpi::recv_coord_data_array[3*i+1];
		atoms::z_coord_array[vmpi::recv_atom_translation_array[i]] = vmpi::recv_coord_data_array[3*i+2];

	}

   #endif

	return;

}

//-----------------------------------------------------------------------------
// Reverse halo exchange for force contributions
//-----------------------------------------------------------------------------
void mpi_sum_halo_forces(std::vector<double>& force_x,
                         std::vector<double>& force_y,
                         std::vector<double>& force_z){
#ifdef MPICF

   // recv_atom_translation_array maps entries received by a normal halo swap onto this rank's halo indices 
   // pack x/y/z from every halo copy and send each group back to its owner
   std::vector<double> send_force_data(3 * vmpi::recv_atom_translation_array.size(), 0.0);
   for(unsigned int i = 0; i < vmpi::recv_atom_translation_array.size(); i++){
      const int atom = vmpi::recv_atom_translation_array[i];
      send_force_data[3*i+0] = force_x[atom];
      send_force_data[3*i+1] = force_y[atom];
      send_force_data[3*i+2] = force_z[atom];
   }

   // send_atom_translation_array maps the owned boundary atoms exported by a normal halo swap
   std::vector<double> recv_force_data(3 * vmpi::send_atom_translation_array.size(), 0.0);
   std::vector<MPI_Request> force_requests;
   force_requests.reserve(2 * vmpi::num_processors);
   MPI_Request req = MPI_REQUEST_NULL;

   for(int p = 0; p < vmpi::num_processors; p++){
      // send forces on local halo copies to their owner, rank p
      if(vmpi::recv_num_array[p] != 0){
         const int num_values = 3 * vmpi::recv_num_array[p];
         const int start = 3 * vmpi::recv_start_index_array[p];
         force_requests.push_back(req);
         MPI_Isend(&send_force_data[start], num_values, MPI_DOUBLE, p, 49,
                   MPI_COMM_WORLD, &force_requests.back());
      }
      // receive forces calculated on rank p's copies of our owned atoms
      if(vmpi::send_num_array[p] != 0){
         const int num_values = 3 * vmpi::send_num_array[p];
         const int start = 3 * vmpi::send_start_index_array[p];
         force_requests.push_back(req);
         MPI_Irecv(&recv_force_data[start], num_values, MPI_DOUBLE, p, 49,
                   MPI_COMM_WORLD, &force_requests.back());
      }
   }

   // wait for all sends and receives to finish, then unpack the received forces
   if(!force_requests.empty()){
      vmpi::TotalComputeTime += vmpi::SwapTimer(vmpi::ComputeTime, vmpi::WaitTime);
      std::vector<MPI_Status> force_status(force_requests.size());
      MPI_Waitall(force_requests.size(), &force_requests[0], &force_status[0]);
      vmpi::TotalWaitTime += vmpi::SwapTimer(vmpi::WaitTime, vmpi::ComputeTime);
   }

   // sum the received forces back onto the owned atoms
   for(unsigned int i = 0; i < vmpi::send_atom_translation_array.size(); i++){
      const int atom = vmpi::send_atom_translation_array[i];
      force_x[atom] += recv_force_data[3*i+0];
      force_y[atom] += recv_force_data[3*i+1];
      force_z[atom] += recv_force_data[3*i+2];
   }
#endif

   return;
}

} // end of namespace vmpi
