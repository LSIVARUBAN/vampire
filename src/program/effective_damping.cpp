//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2014. All rights reserved.
//
//-----------------------------------------------------------------------------

// C++ standard library headers
#include <iostream>
#include <cmath>

// Vampire Header files
#include "atoms.hpp"
#include "errors.hpp"
#include "material.hpp"
#include "program.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "vmath.hpp"
#include "vmpi.hpp"

/// Function to rotate all spins around an axis
void rotate_spins_around_axis(std::vector<double>& axis, double angle_deg){

   double phi = angle_deg * 3.14159265359 / 180.0;
   double u = axis[0];
   double v = axis[1];
   double w = axis[2];
   
   // Normalise axis
   double norm = sqrt(u*u + v*v + w*w);
   if(norm > 0.0){
       u /= norm;
       v /= norm;
       w /= norm;
   }

   double cos_phi = cos(phi);
   double sin_phi = sin(phi);
   double one_minus_cos_phi = 1.0 - cos_phi;

   // loop over all spins and rotate by phi around the rotational axis
   for(int atom = 0; atom < atoms::num_atoms; atom++){

      // Load spin coordinates
      double x = atoms::x_spin_array[atom];
      double y = atoms::y_spin_array[atom];
      double z = atoms::z_spin_array[atom];

      double dot = x*u + y*v + z*w;

      // Calculate rotated spin positions
      double new_x = u * dot * one_minus_cos_phi + x * cos_phi + (-w * y + v * z) * sin_phi;
      double new_y = v * dot * one_minus_cos_phi + y * cos_phi + ( w * x - u * z) * sin_phi;
      double new_z = w * dot * one_minus_cos_phi + z * cos_phi + (-v * x + u * y) * sin_phi;

      atoms::x_spin_array[atom] = new_x;
      atoms::y_spin_array[atom] = new_y;
      atoms::z_spin_array[atom] = new_z;
   }

   return;
}

/// Function to rotate all spin around the x-axis
void rotate_spins_around_x_axis(double ddx){

   std::vector< std::vector<double> > x_rotation_matrix,y_rotation_matrix,z_rotation_matrix;

   // determine rotational matrices for phi, theta rotation
   vmath::set_rotational_matrix(ddx, 0.0, 0.0, x_rotation_matrix,y_rotation_matrix,z_rotation_matrix);

   // Vectors to hold spins
   std::vector<double> Sold(3), Snew(3);

   // loop over all spins and rotate by phi around x
   for(int atom =0;atom<atoms::num_atoms;atom++){

      // Load spin coordinates
      Sold[0]=atoms::x_spin_array[atom];
      Sold[1]=atoms::y_spin_array[atom];
      Sold[2]=atoms::z_spin_array[atom];

      // Calculate new spin positions
      Snew = vmath::matmul(Sold,x_rotation_matrix);

      // Set new spin positions
      atoms::x_spin_array[atom]=Snew[0];
      atoms::y_spin_array[atom]=Snew[1];
      atoms::z_spin_array[atom]=Snew[2];

   }

   return;

}

namespace program{

//-----------------------------------------------------------------------------
//
//   Program to calculate effective gilbert damping from the atomistic model.
//
//   System is first equilibrated at sim:temperature along the z-axis. All
//   spins are then rotated 30 degrees from the z-axis and the magnetisation
//   then relaxes towards the z-direction under the application of a field
//   or magnetocrystalline anisotropy. The magnetization trace can then be
//   fitted to a macroscopic equation to extract the effective damping, alpha.
//
//   Ref. M O A Ellis et al, Physical Review B 86, 174418 (2012)
//
//-----------------------------------------------------------------------------
void effective_damping(){

   // check calling of routine if error checking is activated
   if(err::check==true) std::cout << "program::effective_damping has been called" << std::endl;

   // save system temperature
   double temp=sim::temperature;

   // Set equilibration temperature and field
   sim::temperature=sim::Teq;

   // Equilibrate system
   while(sim::time<sim::equilibration_time){

      sim::integrate(sim::partial_time);

      // Calculate magnetisation statistics
      stats::update();

      // Output data
      vout::data();
   }

   // Rotate all spins around effective_damping_rotation_axis by effective_damping_angle, default is 30 degrees about x-axis
   rotate_spins_around_axis(program::effective_damping_rotation_axis,program::effective_damping_angle);

   // reset system temperature
   sim::temperature=temp;

   // Perform Time Series with relaxation
   while(sim::time<sim::equilibration_time+sim::total_time){

      // Integrate system
      sim::integrate(sim::partial_time);

      // Calculate magnetisation statistics
      stats::update();

      // Output data
      vout::data();

   }

}

}//end of namespace program
