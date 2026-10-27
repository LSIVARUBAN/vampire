//------------------------------------------------------------------------------
//
//   This file is part of the VAMPIRE open source package under the
//   Free BSD licence (see licence file for details).
//
//   (c) Mara Strungaru 2022. All rights reserved.
//
//   Email: mara.strungaru@york.ac.uk
//
//   implementation based on the paper Phys. Rev. B 103, 024429, (2021) M.Strungaru, M.O.A. Ellis et al
//------------------------------------------------------------------------------
//

// C++ standard library headers
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

// Vampire headers
#include "sld.hpp"
#include "atoms.hpp"
#include "create.hpp"
#include "sim.hpp"
#include "errors.hpp"
#include "material.hpp"
#include "vmpi.hpp"
#include "vio.hpp"

// sld module headers
#include "internal.hpp"


namespace sld{

   namespace{

   // if requested output a force table for the first force evaluation containing contributions just from the mechanical potential
   void write_force_debug_table(const int start_index,
                                const int end_index,
                                const std::vector<double>& x_coord_array,
                                const std::vector<double>& y_coord_array,
                                const std::vector<double>& z_coord_array,
                                const std::vector<double>& initial_forces_x,
                                const std::vector<double>& initial_forces_y,
                                const std::vector<double>& initial_forces_z,
                                const std::vector<double>& forces_array_x,
                                const std::vector<double>& forces_array_y,
                                const std::vector<double>& forces_array_z,
                                const std::vector<double>& potential_eng){

      std::ofstream output("force-debug.data");
      if(!output){
         err::zexit("Unable to open force-debug.data");
      }

      output << "# VAMPIRE mechanical-potential force debug from the first force evaluation\n";
      output << "# atom x_A y_A z_A fx_eV_per_A fy_eV_per_A fz_eV_per_A force_magnitude_eV_per_A energy_eV\n";
      output << std::setprecision(17);

      for(int atom = start_index; atom < end_index; atom++){
         const int local_atom = atom - start_index;
         const double fx = forces_array_x[atom] - initial_forces_x[local_atom];
         const double fy = forces_array_y[atom] - initial_forces_y[local_atom];
         const double fz = forces_array_z[atom] - initial_forces_z[local_atom];
         const double magnitude = std::sqrt(fx*fx + fy*fy + fz*fz);

         output << atom << " " << x_coord_array[atom] << " " << y_coord_array[atom] << " " << z_coord_array[atom] << " " << fx << " " << fy << " " << fz << " " << magnitude << " " << potential_eng[atom] << "\n";
      }
   }

   // Phonon pump state for the current velocity half-kick
   std::vector<double> pump_force_x; // pump force on each owned atom (eV/A)
   std::vector<double> pump_force_y;
   std::vector<double> pump_power_before; // F.v for each owned atom before the kick
   bool pump_kick_active = false;
   const double two_pi = 2.0*M_PI;

   // Number of physical atoms owned by this rank, excluding halo copies
   int pump_owned_atoms(){
      #ifdef MPICF
         return vmpi::num_core_atoms + vmpi::num_bdry_atoms;
      #else
         return atoms::num_atoms;
      #endif
   }

   // Time since the start of production (s) at a stage step, as the pulse times are measured from the end of equilibration
   double pump_production_time(const double step){
      return (step - double(sim::equilibration_time))*mp::dt_SI;
   }

   // Half-open pulse window [start,end) in production time
   bool pump_active(const double time){
      const double eps = 4.0*std::numeric_limits<double>::epsilon();
      const double start = internal::phonon_pulse_start_time;
      const double end = internal::phonon_pulse_end_time;
      return time + eps*std::max(std::abs(time),std::abs(start)) >= start && time < end - eps*std::max(std::abs(time),std::abs(end));
   }

   // Travelling wave phase theta = 2 pi f (t - t_start) - q.r + phase, with live positions in Angstrom and the carrier time
   // measured from the start of the pulse, so every frequency starts at the same phase
   double pump_phase(const int atom, const double time){
      const double qr = internal::phonon_wavevector[0]*atoms::x_coord_array[atom] + internal::phonon_wavevector[1]*atoms::y_coord_array[atom] + internal::phonon_wavevector[2]*atoms::z_coord_array[atom];
      return two_pi*internal::phonon_frequency*(time - internal::phonon_pulse_start_time) - qr + internal::phonon_phase;
   }

   // Global means of cos(theta) and sin(theta) over all owned atoms, removed so that the pump exerts no net force
   void pump_mean(const double time, double& mean_cos, double& mean_sin){
      const int num_atoms = pump_owned_atoms();
      double sums[3] = {0.0, 0.0, double(num_atoms)};
      for(int atom=0; atom<num_atoms; ++atom){
         const double theta = pump_phase(atom,time);
         sums[0] += std::cos(theta);
         sums[1] += std::sin(theta);
      }
      #ifdef MPICF
         MPI_Allreduce(MPI_IN_PLACE,sums,3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      #endif
      mean_cos = sums[0]/sums[2];
      mean_sin = sums[1]/sums[2];
   }

   // Pump force (Ax cos(theta), Ay sin(theta), 0) on an atom, minus the global mean
   void pump_force(const int atom, const double time, const double mean_cos, const double mean_sin, double& fx, double& fy){
      const double theta = pump_phase(atom,time);
      fx = internal::phonon_force_amplitude[0]*(std::cos(theta) - mean_cos);
      fy = internal::phonon_force_amplitude[1]*(std::sin(theta) - mean_sin);
   }

   } // end of anonymous namespace

   namespace internal{

      // Phonon pump originally written by Muhammad Hamza Asim
      // Resolve the pump wavevector and check the pump settings
      void initialise_phonon_pump(){
         if(!linear_pump_enabled) return;
         if(sim::integrator != sim::suzuki_trotter) err::zexit("The phonon pump requires sim:integrator = spin-lattice");
         if(phonon_frequency <= 0.0 || phonon_pulse_start_time < 0.0 || phonon_pulse_end_time <= phonon_pulse_start_time) err::zexit("The phonon pump requires phonon:frequency > 0 and 0 <= phonon:pulse-start-time < phonon:pulse-end-time");

         zlog << zTs() << "Phonon pump: frequency " << phonon_frequency << " Hz, pulse " << phonon_pulse_start_time << " to " << phonon_pulse_end_time << " s after the start of production, phase " << phonon_phase << " rad, amplitude (" << phonon_force_amplitude[0] << ", " << phonon_force_amplitude[1] << ") eV/A" << std::endl;

         bool zero_q = true;
         for(int d=0; d<3; ++d){
            // q_d = 2 pi direction_d / lambda_d in 1/A
            phonon_wavevector[d] = phonon_wave_lambda[d] > 0.0 ? two_pi*phonon_wave_direction[d]/phonon_wave_lambda[d] : 0.0;
            // number of wavelengths across the box, which must be an integer along a periodic direction
            const double mode = phonon_wavevector[d]*cs::system_dimensions[d]/two_pi;
            if(cs::pbc[d] && std::abs(mode-std::round(mode)) > 1e-8) err::zexit("The phonon pump wavevector is not suitable with the periodic box");
            zlog << zTs() << "Phonon pump q_" << "xyz"[d] << " = " << phonon_wavevector[d] << " 1/A, box " << cs::system_dimensions[d] << " A, mode index " << mode << std::endl;
            if(phonon_wavevector[d] != 0.0) zero_q = false;
         }
         // a uniform q=0 drive is removed entirely by the mean subtraction
         if(zero_q) err::zexit("The phonon pump requires a nonzero wavevector");
      }

      // Add the pump force to owned atoms just before a velocity half-kick at the given stage step (sim::time or sim::time+1)
      void apply_phonon_pump(const double step){
         const double time = pump_production_time(step);
         pump_kick_active = linear_pump_enabled && pump_active(time);
         if(!pump_kick_active) return;
         const int num_atoms = pump_owned_atoms();
         double mean_cos, mean_sin;
         pump_mean(time,mean_cos,mean_sin);
         pump_force_x.resize(num_atoms);
         pump_force_y.resize(num_atoms);
         pump_power_before.resize(num_atoms);
         for(int atom=0; atom<num_atoms; ++atom){
            double fx, fy;
            pump_force(atom,time,mean_cos,mean_sin,fx,fy);
            pump_force_x[atom] = fx;
            pump_force_y[atom] = fy;
            forces_array_x[atom] += fx;
            forces_array_y[atom] += fy;
            pump_power_before[atom] = fx*atoms::x_velo_array[atom] + fy*atoms::y_velo_array[atom];
         }
      }

      // Accumulate the pump work of the completed half-kick, (dt/2) F.(v_before + v_after)/2, with F in eV/A and v in A/ps
      void finish_phonon_pump_kick(){
         if(!pump_kick_active) return;
         double work = 0.0;
         for(int atom=0; atom<pump_owned_atoms(); ++atom){
            work += pump_power_before[atom] + pump_force_x[atom]*atoms::x_velo_array[atom] + pump_force_y[atom]*atoms::y_velo_array[atom];
         }
         phonon_pump_work += 0.25*mp::dt_SI*1e12*work;
         pump_kick_active = false;
      }

   } // end of internal namespace

   // Update pump work and mode projections at the statistics rate
   void update_phonon_pump_diagnostics(){
      if(!phonon_pump_diagnostics || !internal::linear_pump_enabled) return;
      const int num_atoms = pump_owned_atoms();
      const double ax = internal::phonon_force_amplitude[0];
      const double ay = internal::phonon_force_amplitude[1];
      const double norm = std::sqrt(ax*ax + ay*ay);
      // projections are along the unit vector e = (Ax, Ay, 0)/|A|
      const double inv_norm = norm > 0.0 ? 1.0/norm : 0.0;

      // work, Re/Im of u.e and v.e projected onto exp(-i q.R0), maximum displacement and atom count
      double values[7] = {internal::phonon_pump_work, 0.0, 0.0, 0.0, 0.0, 0.0, double(num_atoms)};
      for(int atom=0; atom<num_atoms; ++atom){
         // displacement from the creation lattice site R0
         const double dx = atoms::x_coord_array[atom] - internal::x0_coord_array[atom];
         const double dy = atoms::y_coord_array[atom] - internal::y0_coord_array[atom];
         const double dz = atoms::z_coord_array[atom] - internal::z0_coord_array[atom];
         const double au = (ax*dx + ay*dy)*inv_norm;
         const double av = (ax*atoms::x_velo_array[atom] + ay*atoms::y_velo_array[atom])*inv_norm;
         const double qr = internal::phonon_wavevector[0]*internal::x0_coord_array[atom] + internal::phonon_wavevector[1]*internal::y0_coord_array[atom] + internal::phonon_wavevector[2]*internal::z0_coord_array[atom];
         const double c = std::cos(qr);
         const double s = -std::sin(qr);
         values[1] += au*c;
         values[2] += au*s;
         values[3] += av*c;
         values[4] += av*s;
         values[5] = std::max(values[5], std::sqrt(dx*dx + dy*dy + dz*dz));
      }
      #ifdef MPICF
         const double max_displacement = values[5];
         MPI_Allreduce(MPI_IN_PLACE,values,7,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
         MPI_Allreduce(&max_displacement,&values[5],1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
      #endif
      for(int i=0; i<6; ++i) phonon_pump_statistics[i] = i == 5 ? values[i] : values[i]/values[6];
   }

   // Append the accumulated pump work as an optional checkpoint record
   void save_phonon_pump_checkpoint(std::ostream& stream){
      if(!internal::linear_pump_enabled) return;
      const uint64_t tag = 0x3230504D55504E50; // PNPUMP02
      stream.write(reinterpret_cast<const char*>(&tag),sizeof(tag));
      stream.write(reinterpret_cast<const char*>(&internal::phonon_pump_work),sizeof(double));
   }

   // Restore the accumulated pump work for a continued run
   void load_phonon_pump_checkpoint(std::istream& stream){
      stream.read(reinterpret_cast<char*>(&internal::phonon_pump_work),sizeof(double));
      // Mechanical potential energies are not in the base checkpoint. Refresh them for the first continued energy sample.
      compute_forces(0,pump_owned_atoms(),atoms::neighbour_list_start_index,atoms::neighbour_list_end_index,
                     atoms::type_array,atoms::neighbour_list_array,internal::x0_coord_array,internal::y0_coord_array,internal::z0_coord_array,
                     atoms::x_coord_array,atoms::y_coord_array,atoms::z_coord_array,
                     internal::forces_array_x,internal::forces_array_y,internal::forces_array_z,internal::potential_eng);
   }

   // A continued pump needs its saved work, otherwise the accumulated work would restart from zero
   void check_phonon_pump_checkpoint(const bool loaded){
      if(internal::linear_pump_enabled && !loaded) err::zexit("Cannot continue a phonon pump without its checkpoint");
   }

   void compute_forces(const int start_index, // first atom for exchange interactions to be calculated
               const int end_index,
               const std::vector<int>& neighbour_list_start_index,
               const std::vector<int>& neighbour_list_end_index,
               const std::vector<int>& type_array, // type for atom
               const std::vector<int>& neighbour_list_array, // list of interactions between atom
               const std::vector<double>& x0_coord_array, // coord vectors for atoms
               const std::vector<double>& y0_coord_array,
               const std::vector<double>& z0_coord_array,
               const std::vector<double>& x_coord_array, // coord vectors for atoms
               const std::vector<double>& y_coord_array,
               const std::vector<double>& z_coord_array,
               std::vector<double>& forces_array_x, //  vectors for forces
               std::vector<double>& forces_array_y,
               std::vector<double>& forces_array_z,
               std::vector<double>& potential_eng){
      const bool write_force_debug = sld::output_force_debug && !sld::internal::force_debug_written;
      std::vector<double> initial_forces_x;
      std::vector<double> initial_forces_y;
      std::vector<double> initial_forces_z;
      // save the forces before the potential adds its contribution, we later subtract these from the final forces to isolate the mechanical potential contribution
      if(write_force_debug){
         const int natoms = end_index - start_index;
         initial_forces_x.resize(natoms);
         initial_forces_y.resize(natoms);
         initial_forces_z.resize(natoms);
         for(int atom = start_index; atom < end_index; atom++){
            const int local_atom = atom - start_index;
            initial_forces_x[local_atom] = forces_array_x[atom];
            initial_forces_y[local_atom] = forces_array_y[atom];
            initial_forces_z[local_atom] = forces_array_z[atom];
         }
      }

      switch(sld::internal::lattice_potential){
         case sld::internal::harmonic_lattice_potential:
            // Harmonic and Morse are simple pair potentials
            // the harmonic reference distances are precomputed from the x0, y0 and z0 coordinates at initialisation
            internal::compute_forces_harmonic(start_index, end_index,
                                              neighbour_list_start_index, neighbour_list_end_index,
                                              type_array, neighbour_list_array,
                                              sld::internal::harmonic_reference_distance_array,
                                              x_coord_array, y_coord_array, z_coord_array,
                                              forces_array_x, forces_array_y, forces_array_z, potential_eng);
            break;

         case sld::internal::morse_lattice_potential:
            internal::compute_forces_morse(start_index, end_index,
                                           neighbour_list_start_index, neighbour_list_end_index,
                                           type_array, neighbour_list_array,
                                           x_coord_array, y_coord_array, z_coord_array,
                                           forces_array_x, forces_array_y, forces_array_z, potential_eng);
            break;

         case sld::internal::snap_lattice_potential:
            // SNAP is many body therefore the force on atom i depends on the neighbour density around each central atom (not just a pair)
            internal::compute_forces_snap(start_index, end_index,
                                          neighbour_list_start_index, neighbour_list_end_index,
                                          type_array, neighbour_list_array,
                                          x_coord_array, y_coord_array, z_coord_array,
                                          forces_array_x, forces_array_y, forces_array_z, potential_eng);
            break;

         case sld::internal::snap_zbl_lattice_potential:
            // SNAP plus a ZBL potential overlay for short range pair repulsion
            internal::compute_forces_snap_zbl(start_index, end_index,
                                             neighbour_list_start_index, neighbour_list_end_index,
                                             type_array, neighbour_list_array,
                                             x_coord_array, y_coord_array, z_coord_array,
                                             forces_array_x, forces_array_y, forces_array_z, potential_eng);
            break;

         default:
            break;
      }

      // write forces from mechnical potential to file once, if requested
      if(write_force_debug){
         write_force_debug_table(start_index, end_index,
                                 x_coord_array, y_coord_array, z_coord_array,
                                 initial_forces_x, initial_forces_y, initial_forces_z,
                                 forces_array_x, forces_array_y, forces_array_z,
                                 potential_eng);
         sld::internal::force_debug_written = true;
      }
      
      return;

}

namespace internal{

namespace{

struct zbl_coeff_t{
   double d1a; // first exponential decay constant divided by screening length
   double d2a; // second
   double d3a; // third
   double d4a; // fourth
   double zze; // Zi*Zj*e^2/(4*pi*epsilon0) in eV Angstrom
   double sw1; // force switch coefficient A
   double sw2; // force switch coefficient B
   double sw3; // energy switch coefficient A/3
   double sw4; // energy switch coefficient B/4
   double sw5; // energy constant offset
};

double e_zbl(const double r, const zbl_coeff_t& coeff){

   // unswitched screened Coulomb pair energy in eV
   const double c1 = 0.02817;
   const double c2 = 0.28022;
   const double c3 = 0.50986;
   const double c4 = 0.18175;
   const double d1aij = coeff.d1a;
   const double d2aij = coeff.d2a;
   const double d3aij = coeff.d3a;
   const double d4aij = coeff.d4a;
   const double zzeij = coeff.zze;
   const double rinv = 1.0 / r;

   // phi(r/a) is a weighted sum of four exponentials with each term describing one part of the ZBL screening cloud
   double sum = c1 * std::exp(-d1aij * r);
   sum += c2 * std::exp(-d2aij * r);
   sum += c3 * std::exp(-d3aij * r);
   sum += c4 * std::exp(-d4aij * r);

   // E_ZBL(r) = Zi*Zj*e^2/(4*pi*epsilon0) * phi(r/a) / r (in eV)
   double result = zzeij * sum * rinv;

   return result;

}

double dzbldr(const double r, const zbl_coeff_t& coeff){

   // dE_ZBL/dr before the switching derivative is added
   const double c1 = 0.02817;
   const double c2 = 0.28022;
   const double c3 = 0.50986;
   const double c4 = 0.18175;
   const double d1aij = coeff.d1a;
   const double d2aij = coeff.d2a;
   const double d3aij = coeff.d3a;
   const double d4aij = coeff.d4a;
   const double zzeij = coeff.zze;
   const double rinv = 1.0 / r;

   const double e1 = std::exp(-d1aij * r);
   const double e2 = std::exp(-d2aij * r);
   const double e3 = std::exp(-d3aij * r);
   const double e4 = std::exp(-d4aij * r);

   // phi(r/a)
   double sum = c1 * e1;
   sum += c2 * e2;
   sum += c3 * e3;
   sum += c4 * e4;

   // d(phi)/dr
   double sum_p = -c1 * d1aij * e1;
   sum_p -= c2 * d2aij * e2;
   sum_p -= c3 * d3aij * e3;
   sum_p -= c4 * d4aij * e4;

   // dE/dr = zze * (phi'/r - phi/r^2).
   double result = zzeij * (sum_p - sum * rinv) * rinv;

   return result;

}

double d2zbldr2(const double r, const zbl_coeff_t& coeff){

   const double c1 = 0.02817;
   const double c2 = 0.28022;
   const double c3 = 0.50986;
   const double c4 = 0.18175;
   const double d1aij = coeff.d1a;
   const double d2aij = coeff.d2a;
   const double d3aij = coeff.d3a;
   const double d4aij = coeff.d4a;
   const double zzeij = coeff.zze;
   const double rinv = 1.0 / r;

   const double e1 = std::exp(-d1aij * r);
   const double e2 = std::exp(-d2aij * r);
   const double e3 = std::exp(-d3aij * r);
   const double e4 = std::exp(-d4aij * r);

   // phi(r/a)
   double sum = c1 * e1;
   sum += c2 * e2;
   sum += c3 * e3;
   sum += c4 * e4;

   // -d(phi)/dr
   double sum_p = c1 * e1 * d1aij;
   sum_p += c2 * e2 * d2aij;
   sum_p += c3 * e3 * d3aij;
   sum_p += c4 * e4 * d4aij;

   // d2(phi)/dr2
   double sum_pp = c1 * e1 * d1aij * d1aij;
   sum_pp += c2 * e2 * d2aij * d2aij;
   sum_pp += c3 * e3 * d3aij * d3aij;
   sum_pp += c4 * e4 * d4aij * d4aij;

   // d2E/dr2 for E = zze * phi/r.
   double result = zzeij * (sum_pp + 2.0 * sum_p * rinv + 2.0 * sum * rinv * rinv) * rinv;

   return result;

}

zbl_coeff_t set_zbl_coeff(const double zi,
                          const double zj,
                          const double cut_inner,
                          const double cut_outer){

   const double pzbl = 0.23;
   const double a0 = 0.46850;
   const double d1 = 0.20162;
   const double d2 = 0.40290;
   const double d3 = 0.94229;
   const double d4 = 3.19980;
   const double qqr2e_metal = 14.399645;

   zbl_coeff_t coeff;
   // ZBL screening length: a = a0 / (Zi^0.23 + Zj^0.23)
   const double ainv = (std::pow(zi, pzbl) + std::pow(zj, pzbl)) / a0;
   coeff.d1a = d1 * ainv;
   coeff.d2a = d2 * ainv;
   coeff.d3a = d3 * ainv;
   coeff.d4a = d4 * ainv;
   // Coulomb prefactor
   coeff.zze = zi * zj * qqr2e_metal;

   // e      = t^3 * (sw3 + sw4*t) + sw5
   // dedr   = t^2 * (sw1 + sw2*t)
   // d2edr2 = 2*sw1*t + 3*sw2*t^2
   const double tc = cut_outer - cut_inner;
   const double fc = e_zbl(cut_outer, coeff);
   const double fcp = dzbldr(cut_outer, coeff);
   const double fcpp = d2zbldr2(cut_outer, coeff);

   const double swa = (-3.0 * fcp + tc * fcpp) / (tc * tc);
   const double swb = (2.0 * fcp - tc * fcpp) / (tc * tc * tc);
   const double swc = -fc + (tc / 2.0) * fcp - (tc * tc / 12.0) * fcpp;

   coeff.sw1 = swa;
   coeff.sw2 = swb;
   coeff.sw3 = swa / 3.0;
   coeff.sw4 = swb / 4.0;
   coeff.sw5 = swc;

   return coeff;

}

double zbl_single_energy(const double r,
                         const double cut_inner,
                         const zbl_coeff_t& coeff){

   double phi = e_zbl(r, coeff);
   phi += coeff.sw5;

   if(r > cut_inner){
      const double t = r - cut_inner;
      const double eswitch = t*t*t * (coeff.sw3 + coeff.sw4 * t);
      phi += eswitch;
   }

   return phi;

}

// return dE/dr (not force)
double zbl_single_dedr(const double r,
                       const double cut_inner,
                       const zbl_coeff_t& coeff){

   double fforce = dzbldr(r, coeff);

   if(r > cut_inner){
      const double t = r - cut_inner;
      const double fswitch = t*t * (coeff.sw1 + coeff.sw2 * t);
      fforce += fswitch;
   }

   return fforce;

}

std::vector<int> count_neighbours_within_cutoff(const int start_index,
                                                const int end_index,
                                                const std::vector<int>& neighbour_list_start_index,
                                                const std::vector<int>& neighbour_list_end_index,
                                                const std::vector<int>& neighbour_list_array,
                                                const std::vector<double>& x_coord_array,
                                                const std::vector<double>& y_coord_array,
                                                const std::vector<double>& z_coord_array,
                                                const double cutoff){

   std::vector<int> neighbour_counts(end_index - start_index, 0);
   const double cutoff_sqr = cutoff * cutoff;

   for(int i = start_index; i < end_index; i++){
      const double rx = x_coord_array[i];
      const double ry = y_coord_array[i];
      const double rz = z_coord_array[i];
      const int nbr_start = neighbour_list_start_index[i];
      const int nbr_end = neighbour_list_end_index[i] + 1;

      for(int n = nbr_start; n < nbr_end; n++){
         const int j = neighbour_list_array[n];
         if(j == i) continue;

         double dx = rx - x_coord_array[j];
         double dy = ry - y_coord_array[j];
         double dz = rz - z_coord_array[j];

         dx = sld::PBC_wrap(dx, cs::system_dimensions[0], cs::pbc[0]);
         dy = sld::PBC_wrap(dy, cs::system_dimensions[1], cs::pbc[1]);
         dz = sld::PBC_wrap(dz, cs::system_dimensions[2], cs::pbc[2]);

         const double rsq = dx*dx + dy*dy + dz*dz;
         if(rsq > 0.0 && rsq < cutoff_sqr) neighbour_counts[i - start_index]++;
      }
   }

   return neighbour_counts;

}

} // end of anonymous namespace

void print_force_debug_summary(const std::string& label,
                               const int call,
                               const int start_index,
                               const int end_index,
                               const std::vector<int>& neighbour_counts,
                               const std::vector<double>& initial_forces_x,
                               const std::vector<double>& initial_forces_y,
                               const std::vector<double>& initial_forces_z,
                               const std::vector<double>& forces_array_x,
                               const std::vector<double>& forces_array_y,
                               const std::vector<double>& forces_array_z){

   const int natoms = end_index - start_index;
   if(natoms <= 0) return;

   int min_neigh = std::numeric_limits<int>::max();
   int max_neigh = 0;
   double avg_neigh = 0.0;
   double min_force = std::numeric_limits<double>::max();
   double max_force = 0.0;
   double avg_force = 0.0;
   int nonfinite_force = 0;

   for(int atom = start_index; atom < end_index; atom++){
      const int local_atom = atom - start_index;
      const int ninside = neighbour_counts[local_atom];
      min_neigh = std::min(min_neigh, ninside);
      max_neigh = std::max(max_neigh, ninside);
      avg_neigh += double(ninside);

      const double fx = forces_array_x[atom] - initial_forces_x[local_atom];
      const double fy = forces_array_y[atom] - initial_forces_y[local_atom];
      const double fz = forces_array_z[atom] - initial_forces_z[local_atom];
      const double force = std::sqrt(fx*fx + fy*fy + fz*fz);

      if(std::isfinite(force)){
         min_force = std::min(min_force, force);
         max_force = std::max(max_force, force);
         avg_force += force;
      }
      else nonfinite_force++;
   }

   avg_neigh /= double(natoms);
   const int finite_force_count = natoms - nonfinite_force;
   if(finite_force_count > 0) avg_force /= double(finite_force_count);
   else min_force = 0.0;

   const std::streamsize old_precision = std::cout.precision();
   std::cout<<std::setprecision(12);
   std::cout<<label<<" force debug call "<<call<<std::endl;
   std::cout<<"neighbours avg/min/max [count]: "<<avg_neigh<<" "<<min_neigh<<" "<<max_neigh<<std::endl;
   std::cout<<"force avg/min/max [eV/Angstrom]: "<<avg_force<<" "<<min_force<<" "<<max_force<<std::endl;
   if(nonfinite_force > 0){
      std::cout<<"nonfinite force count: "<<nonfinite_force<<std::endl;
   }

   std::vector<int> sample_atoms;
   sample_atoms.push_back(start_index);
   if(natoms > 2) sample_atoms.push_back(start_index + natoms / 2);
   if(natoms > 1) sample_atoms.push_back(end_index - 1);

   for(size_t isample = 0; isample < sample_atoms.size(); isample++){
      const int atom = sample_atoms[isample];
      const int local_atom = atom - start_index;
      const double fx = forces_array_x[atom] - initial_forces_x[local_atom];
      const double fy = forces_array_y[atom] - initial_forces_y[local_atom];
      const double fz = forces_array_z[atom] - initial_forces_z[local_atom];
      const double force = std::sqrt(fx*fx + fy*fy + fz*fz);

      std::cout<<"atom "<<atom
               <<" neighbours [count]: "<<neighbour_counts[local_atom]
               <<" force [eV/Angstrom]: "
               <<fx<<" "<<fy<<" "<<fz
               <<" magnitude "<<force<<std::endl;
   }

   std::cout.precision(old_precision);

}



// Store the minimum image reference distance r_ij0 [A] for every neighbour list entry at initiisation since x0, y0 and z0 are fixed (avoids a pbc wrap for every harmonic force evaluation)
void initialise_harmonic_reference_distances(){

   harmonic_reference_distance_array.assign(atoms::neighbour_list_array.size(), 0.0);

   for(int i = 0; i < atoms::num_atoms; i++){
      const int nbr_start = atoms::neighbour_list_start_index[i];
      const int nbr_end = atoms::neighbour_list_end_index[i] + 1;

      for(int n = nbr_start; n < nbr_end; n++){
         const int j = atoms::neighbour_list_array[n];
         if(j == i) continue; // i=j skipped by the force loop so leave as zero

         double dx0 = x0_coord_array[i] - x0_coord_array[j];
         double dy0 = y0_coord_array[i] - y0_coord_array[j];
         double dz0 = z0_coord_array[i] - z0_coord_array[j];

         dx0 = sld::PBC_wrap(dx0, cs::system_dimensions[0], cs::pbc[0]);
         dy0 = sld::PBC_wrap(dy0, cs::system_dimensions[1], cs::pbc[1]);
         dz0 = sld::PBC_wrap(dz0, cs::system_dimensions[2], cs::pbc[2]);

         harmonic_reference_distance_array[n] = sqrt(dx0*dx0 + dy0*dy0 + dz0*dz0);
      }
   }

}

template<bool fully_periodic>
void compute_forces_harmonic(const int start_index,
            const int end_index, // last +1 atom to be calculated
            const std::vector<int>& neighbour_list_start_index,
            const std::vector<int>& neighbour_list_end_index,
            const std::vector<int>& type_array, // type for atom
            const std::vector<int>& neighbour_list_array, // list of interactions between atom
            const std::vector<double>& reference_distance_array, // reference distance for each neighbour list entry
            const std::vector<double>& x_coord_array, // coord vectors for atoms
            const std::vector<double>& y_coord_array,
            const std::vector<double>& z_coord_array,
            std::vector<double>& forces_array_x, //  vectors for forces
            std::vector<double>& forces_array_y,
            std::vector<double>& forces_array_z,
            std::vector<double>& potential_eng){



            double rx, ry, rz;
            double dx, dy, dz;
            double fx = 0.0, fy = 0.0, fz = 0.0;
            double rji_sqr, rji, rji0, inv_rji;
            int j, total_int;
            double r_sqr_cut=sld::internal::r_cut_pot*sld::internal::r_cut_pot;
            double energy;
            const int natoms = end_index - start_index;
            const bool debug_this_call = sld::internal::harmonic_debug_enabled &&
                                         sld::internal::harmonic_debug_force_calls <
                                         sld::internal::harmonic_debug_max_force_calls;
            std::vector<int> neighbour_counts;
            std::vector<double> initial_forces_x;
            std::vector<double> initial_forces_y;
            std::vector<double> initial_forces_z;

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



            for(int i=start_index;i< end_index; ++i){

               fx = 0.0;
               fy = 0.0;
               fz = 0.0;
               energy=0.0;
               total_int=0;

               rx = x_coord_array[i];
               ry = y_coord_array[i];
               rz = z_coord_array[i];

                //note for sld_neighbour_list_array
                // int nbr_end = neighbour_list_end_index[i];
                //for(int i=start_index;i<= end_index; ++i)
                //for( int n = nbr_start; n <=nbr_end; ++n)



        	      int nbr_start = neighbour_list_start_index[i];
        	      int nbr_end = neighbour_list_end_index[i]+1;

        	      for( int n = nbr_start; n < nbr_end; ++n){
        	        j = neighbour_list_array[n];

        	        if ( j != i){
        		       dx = -x_coord_array[j] + rx;
        		       dy = -y_coord_array[j] + ry;
        		       dz = -z_coord_array[j] + rz;

                   dx = sld::PBC_wrap<fully_periodic>(dx, cs::system_dimensions[0], cs::pbc[0]);
                   dy = sld::PBC_wrap<fully_periodic>(dy, cs::system_dimensions[1], cs::pbc[1]);
                   dz = sld::PBC_wrap<fully_periodic>(dz, cs::system_dimensions[2], cs::pbc[2]);

        		       rji_sqr = dx*dx + dy*dy + dz*dz;

        		       if( rji_sqr < r_sqr_cut){

                       total_int++;


        		           rji = sqrt(rji_sqr);
                       rji0 = reference_distance_array[n]; // fixed reference distance stored at initialisation
        		           inv_rji = 1.0/ rji;

                       energy += (rji-rji0)*(rji-rji0);


		               fx -=  (rji-rji0)*dx*inv_rji ; //2 (rji-rj0)*dx*inv_rji -> 2 went at the end
                       fy -=  (rji-rji0)*dy*inv_rji ;
                       fz -=  (rji-rji0)*dz*inv_rji ;


                     }
        	   	}
        	    }
             const unsigned int imat = type_array[i];

             double V0=sld::internal::mp[imat].V0.get(); //0.15

        	    forces_array_x[i] += V0 * 2.0 * fx;
        	    forces_array_y[i] += V0  * 2.0 * fy;
        	    forces_array_z[i] += V0  * 2.0 * fz;
                potential_eng[i] = 0.5 * V0 * energy;

                if(debug_this_call){
                   const int local_atom = i - start_index;
                   neighbour_counts[local_atom] = total_int;
                }

  }

            if(debug_this_call){
               sld::internal::print_force_debug_summary("Harmonic",
                                                        sld::internal::harmonic_debug_force_calls + 1,
                                                        start_index, end_index,
                                                        neighbour_counts,
                                                        initial_forces_x,
                                                        initial_forces_y,
                                                        initial_forces_z,
                                                        forces_array_x,
                                                        forces_array_y,
                                                        forces_array_z);
               sld::internal::harmonic_debug_force_calls++;
            }


     return;
            }

void compute_forces_harmonic(const int start_index,
            const int end_index,
            const std::vector<int>& neighbour_list_start_index,
            const std::vector<int>& neighbour_list_end_index,
            const std::vector<int>& type_array,
            const std::vector<int>& neighbour_list_array,
            const std::vector<double>& reference_distance_array,
            const std::vector<double>& x_coord_array,
            const std::vector<double>& y_coord_array,
            const std::vector<double>& z_coord_array,
            std::vector<double>& forces_array_x,
            std::vector<double>& forces_array_y,
            std::vector<double>& forces_array_z,
            std::vector<double>& potential_eng){

   // Select the PBC mode once to avoid repeated boundary checks in the inner loops
   const bool fully_periodic_boundaries = cs::pbc[0] && cs::pbc[1] && cs::pbc[2];
   if(fully_periodic_boundaries){
      compute_forces_harmonic<true>(
         start_index, end_index, neighbour_list_start_index,
         neighbour_list_end_index, type_array, neighbour_list_array,
         reference_distance_array,
         x_coord_array, y_coord_array, z_coord_array,
         forces_array_x, forces_array_y, forces_array_z, potential_eng);
   }
   else{
      compute_forces_harmonic<false>(
         start_index, end_index, neighbour_list_start_index,
         neighbour_list_end_index, type_array, neighbour_list_array,
         reference_distance_array,
         x_coord_array, y_coord_array, z_coord_array,
         forces_array_x, forces_array_y, forces_array_z, potential_eng);
   }
}

template<bool fully_periodic>
void compute_forces_morse(const int start_index,
            const int end_index, // last +1 atom to be calculated
            const std::vector<int>& neighbour_list_start_index,
            const std::vector<int>& neighbour_list_end_index,
            const std::vector<int>& type_array, // type for atom
            const std::vector<int>& neighbour_list_array, // list of interactions between atom
            const std::vector<double>& x_coord_array, // coord vectors for atoms
            const std::vector<double>& y_coord_array,
            const std::vector<double>& z_coord_array,
            std::vector<double>& forces_array_x, //  vectors for forces
            std::vector<double>& forces_array_y,
            std::vector<double>& forces_array_z,
            std::vector<double>& potential_eng){



            double rx, ry, rz;
            double dx, dy, dz;
            double fx = 0.0, fy = 0.0, fz = 0.0;
            double rji_sqr, rji, inv_rji; // rji0,
            int j; //, total_int;
            double r_sqr_cut=sld::internal::r_cut_pot*sld::internal::r_cut_pot;
            double energy;

           double alpha_m= sld::internal::alpha_m; //1.3885;
           //double r0_m= sld::internal::r0_m;//2.845;
           double morse_D=sld::internal::morse_D; //0.4174;;
           double morse_beta=sld::internal::morse_beta;//exp( alpha_m * r0_m);
           double morse_factor = sld::internal::morse_factor; //-2.0 * morse_D * alpha_m;




            for(int i=start_index;i< end_index; ++i){

               fx = 0.0;
               fy = 0.0;
               fz = 0.0;
               energy=0.0;
               //total_int=0;


               rx = x_coord_array[i];
               ry = y_coord_array[i];
               rz = z_coord_array[i];


        	      int nbr_start = neighbour_list_start_index[i];
        	      int nbr_end = neighbour_list_end_index[i]+1;

        	      for( int n = nbr_start; n < nbr_end; ++n){
        	        j = neighbour_list_array[n];

        	        if ( j != i){
        		       dx = x_coord_array[j] -rx;
        		       dy = y_coord_array[j]- ry;
        		       dz = z_coord_array[j]- rz;

                   dx = sld::PBC_wrap<fully_periodic>( dx, cs::system_dimensions[0], cs::pbc[0]);
                   dy = sld::PBC_wrap<fully_periodic>( dy, cs::system_dimensions[1], cs::pbc[1]);
                   dz = sld::PBC_wrap<fully_periodic>( dz, cs::system_dimensions[2], cs::pbc[2]);



        		       rji_sqr = dx*dx + dy*dy + dz*dz;

        		       if( rji_sqr < r_sqr_cut){

        		           rji = sqrt(rji_sqr);
    		               inv_rji = 1.0/ rji;

		                   double y = morse_beta * exp( - alpha_m * rji);
  		                   double f_morse = y * ( y - 1.0);

        		           //std::cout<<"embedded "<<i<<"\t"<<j<<"\t"<<rji<<"\t"<<position_rij<<"\t"<<int(position_rij)<<"\t"<<v_rij<<"\t"<<rho_rij<<std::endl;



		                 fx +=  f_morse*dx*inv_rji;
                         fy +=  f_morse*dy*inv_rji;
                         fz +=  f_morse*dz*inv_rji;

                       //if(i==0) std::cout<<"fxyz "<<i <<"\t"<<j<<"\t"<<rji<<"\t"<<f_morse*dx*inv_rji<<"\t"<<f_morse*dy*inv_rji<<"\t"<<f_morse*dz*inv_rji<< std::endl;//<<(rji-rji0)*dx*inv_rji<<"\t"<<(rji-rji0)*dy*inv_rji<<"\t"<<(rji-rji0)*dz*inv_rji<<std::endl;
                       //if(i==0) std::cout<<"fxyz "<<i <<"\t"<<j<<"\t"<<rji<<"\t"<<x_coord_array[j]<<"\t"<<y_coord_array[j]<<"\t"<<z_coord_array[j]<< "\t"<<dx<<"\t"<<dy<<"\t"<<dz<<"\t"<<f_morse<<std::endl;//<<(rji-rji0)*dx*inv_rji<<"\t"<<(rji-rji0)*dy*inv_rji<<"\t"<<(rji-rji0)*dz*inv_rji<<std::endl;

                       //if(i==100) std::cout<<"fxyz "<<i <<"\t"<<j<<"\t"<<type_array[i]<<"\t"<<inv_rji<<std::endl;//<<(rji-rji0)*dx*inv_rji<<"\t"<<(rji-rji0)*dy*inv_rji<<"\t"<<(rji-rji0)*dz*inv_rji<<std::endl;
                         energy += y * ( y - 2.0);

                     }


        	   	}


        	    }





        	    forces_array_x[i] += fx*morse_factor;
        	    forces_array_y[i] += fy*morse_factor;
        	    forces_array_z[i] += fz*morse_factor;

              potential_eng[i] = morse_D * energy;

  }



   return;
          }

void compute_forces_morse(const int start_index,
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
            std::vector<double>& potential_eng){

   // Select the PBC mode once to avoid repeated boundary checks in the inner loops
   const bool fully_periodic_boundaries =
      cs::pbc[0] && cs::pbc[1] && cs::pbc[2];
   if(fully_periodic_boundaries){
      compute_forces_morse<true>(
         start_index, end_index, neighbour_list_start_index,
         neighbour_list_end_index, type_array, neighbour_list_array,
         x_coord_array, y_coord_array, z_coord_array,
         forces_array_x, forces_array_y, forces_array_z, potential_eng);
   }
   else{
      compute_forces_morse<false>(
         start_index, end_index, neighbour_list_start_index,
         neighbour_list_end_index, type_array, neighbour_list_array,
         x_coord_array, y_coord_array, z_coord_array,
         forces_array_x, forces_array_y, forces_array_z, potential_eng);
   }
}


void compute_forces_snap(const int start_index,
            const int end_index, // last +1 atom to be calculated
            const std::vector<int>& neighbour_list_start_index,
            const std::vector<int>& neighbour_list_end_index,
            const std::vector<int>& type_array, // type for atom
            const std::vector<int>& neighbour_list_array, // list of interactions between atom
            const std::vector<double>& x_coord_array, // coord vectors for atoms
            const std::vector<double>& y_coord_array,
            const std::vector<double>& z_coord_array,
            std::vector<double>& forces_array_x, //  vectors for forces
            std::vector<double>& forces_array_y,
            std::vector<double>& forces_array_z,
            std::vector<double>& potential_eng){

   sld::internal::snap_potential.compute_forces(start_index, end_index,
                                                neighbour_list_start_index, neighbour_list_end_index,
                                                type_array, neighbour_list_array,
                                                x_coord_array, y_coord_array, z_coord_array,
                                                forces_array_x, forces_array_y, forces_array_z,
                                                potential_eng);

   return;
}

template<bool fully_periodic>
void compute_forces_zbl_overlay(const int start_index,
            const int end_index, // last +1 atom to be calculated
            const std::vector<int>& neighbour_list_start_index,
            const std::vector<int>& neighbour_list_end_index,
            const std::vector<int>& type_array, // type for atom
            const std::vector<int>& neighbour_list_array, // list of interactions between atom
            const std::vector<double>& x_coord_array, // coord vectors for atoms
            const std::vector<double>& y_coord_array,
            const std::vector<double>& z_coord_array,
            std::vector<double>& forces_array_x, //  vectors for forces
            std::vector<double>& forces_array_y,
            std::vector<double>& forces_array_z,
            std::vector<double>& potential_eng){

   // cut_inner is where the smooth switch begins, cut_outer is where the ZBL contribution has gone to zero
   const double cut_inner = sld::internal::zbl_inner_cutoff;
   const double cut_outer = sld::internal::zbl_outer_cutoff;
   const double cut_outer_sqr = cut_outer * cut_outer;
   const double z = sld::internal::zbl_atomic_number;
   const zbl_coeff_t coeff = set_zbl_coeff(z, z, cut_inner, cut_outer);

   // the spin lattice neighbour list can contain repeated entries or periodic/MPI copies of the same physical neighbour
   // so filter by the physical identity rather than local index to ensure ZBL contributions are not added multiple times for the same neighbour
   std::vector<int> seen_atom_identity(atoms::num_atoms, -1);
   int central_atom_marker = 0;

   // loop over central atoms owned by this force call
   // the arrays contain field or SNAP contributions so accumulate with the zbl forces
   for(int i = start_index; i < end_index; i++){
      const double rx = x_coord_array[i];
      const double ry = y_coord_array[i];
      const double rz = z_coord_array[i];

      // accumulate only the ZBL force contribution for this central atom
      double fx = 0.0;
      double fy = 0.0;
      double fz = 0.0;
      double energy = 0.0;
      const int central_atom_identity = atoms::physical_atom_id_array[i];

      // use the spin-lattice neighbour list
      const int nbr_start = neighbour_list_start_index[i];
      const int nbr_end = neighbour_list_end_index[i] + 1;

      for(int n = nbr_start; n < nbr_end; n++){
         const int j = neighbour_list_array[n];
         const int neighbour_atom_identity = atoms::physical_atom_id_array[j];

         // exclude the central atom's own periodic/halo copies aswell as neighbours already visited for this central atom
         if(neighbour_atom_identity != central_atom_identity && seen_atom_identity[neighbour_atom_identity] != central_atom_marker){

            seen_atom_identity[neighbour_atom_identity] = central_atom_marker; // mark this neighbour as visited for this central atom

            // vector from neighbour j to central atom i
            double dx = rx - x_coord_array[j];
            double dy = ry - y_coord_array[j];
            double dz = rz - z_coord_array[j];

            dx = sld::PBC_wrap<fully_periodic>(dx, cs::system_dimensions[0], cs::pbc[0]);
            dy = sld::PBC_wrap<fully_periodic>(dy, cs::system_dimensions[1], cs::pbc[1]);
            dz = sld::PBC_wrap<fully_periodic>(dz, cs::system_dimensions[2], cs::pbc[2]);

            const double rsq = dx*dx + dy*dy + dz*dz;

            // only pairs inside the ZBL outer cutoff contribute
            if(rsq > 0.0 && rsq < cut_outer_sqr){
               const double r = std::sqrt(rsq);

               // fpair = -(dE/dr)/r, then multiply by delx, dely, delz to get the force
               const double fpair = -zbl_single_dedr(r, cut_inner, coeff) / r;

               // Cartesian force on central atom i from neighbour j
               fx += dx * fpair;
               fy += dy * fpair;
               fz += dz * fpair;

               // the central atom loop visits i->j and j->i, so each pair is counted twice, so half the zbl energy here
               energy += 0.5 * zbl_single_energy(r, cut_inner, coeff);
            }
         }
      }
      central_atom_marker++;

      // overlay the ZBL contribution onto any previous force terms
      forces_array_x[i] += fx;
      forces_array_y[i] += fy;
      forces_array_z[i] += fz;
      potential_eng[i] += energy;
   }

}

void compute_forces_zbl_overlay(const int start_index,
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
            std::vector<double>& potential_eng){

   // ZBL is evaluated for every unique neighbour pair
   // Select the PBC mode once to avoid repeated boundary checks in the inner loops
   const bool fully_periodic_boundaries =
      cs::pbc[0] && cs::pbc[1] && cs::pbc[2];
   if(fully_periodic_boundaries){
      compute_forces_zbl_overlay<true>(
         start_index, end_index, neighbour_list_start_index,
         neighbour_list_end_index, type_array, neighbour_list_array,
         x_coord_array, y_coord_array, z_coord_array,
         forces_array_x, forces_array_y, forces_array_z, potential_eng);
   }
   else{
      compute_forces_zbl_overlay<false>(
         start_index, end_index, neighbour_list_start_index,
         neighbour_list_end_index, type_array, neighbour_list_array,
         x_coord_array, y_coord_array, z_coord_array,
         forces_array_x, forces_array_y, forces_array_z, potential_eng);
   }
}

void compute_forces_snap_zbl(const int start_index,
            const int end_index, // last +1 atom to be calculated
            const std::vector<int>& neighbour_list_start_index,
            const std::vector<int>& neighbour_list_end_index,
            const std::vector<int>& type_array, // type for atom
            const std::vector<int>& neighbour_list_array, // list of interactions between atom
            const std::vector<double>& x_coord_array, // coord vectors for atoms
            const std::vector<double>& y_coord_array,
            const std::vector<double>& z_coord_array,
            std::vector<double>& forces_array_x, //  vectors for forces
            std::vector<double>& forces_array_y,
            std::vector<double>& forces_array_z,
            std::vector<double>& potential_eng){

   const int natoms = end_index - start_index;
   const bool debug_this_call = sld::internal::snap_potential.debug_this_call();
   std::vector<int> neighbour_counts;
   std::vector<double> initial_forces_x;
   std::vector<double> initial_forces_y;
   std::vector<double> initial_forces_z;

   if(debug_this_call){
      neighbour_counts = count_neighbours_within_cutoff(start_index, end_index,
                                                        neighbour_list_start_index,
                                                        neighbour_list_end_index,
                                                        neighbour_list_array,
                                                        x_coord_array,
                                                        y_coord_array,
                                                        z_coord_array,
                                                        sld::internal::r_cut_pot);
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

   // calculate SNAP and fill potential_eng with SNAP atomic energies and add SNAP forces
   sld::internal::snap_potential.compute_forces(start_index, end_index,
                                                neighbour_list_start_index, neighbour_list_end_index,
                                                type_array, neighbour_list_array,
                                                x_coord_array, y_coord_array, z_coord_array,
                                                forces_array_x, forces_array_y, forces_array_z,
                                                potential_eng,
                                                false);

   // add the pairwise ZBL repulsion
   sld::internal::compute_forces_zbl_overlay(start_index, end_index,
                                             neighbour_list_start_index, neighbour_list_end_index,
                                             type_array, neighbour_list_array,
                                             x_coord_array, y_coord_array, z_coord_array,
                                             forces_array_x, forces_array_y, forces_array_z,
                                             potential_eng);

   if(debug_this_call){
      sld::internal::print_force_debug_summary("SNAP+ZBL",
                                               sld::internal::snap_potential.debug_call_number(),
                                               start_index, end_index,
                                               neighbour_counts,
                                               initial_forces_x,
                                               initial_forces_y,
                                               initial_forces_z,
                                               forces_array_x,
                                               forces_array_y,
                                               forces_array_z);
      sld::internal::snap_potential.increment_debug_force_calls();
   }

}

         } //end of internal
      } // end of sld namespace
