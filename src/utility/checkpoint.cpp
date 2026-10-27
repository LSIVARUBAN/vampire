//-----------------------------------------------------------------------------
//
// This source file is part of the VAMPIRE open source package under the
// GNU GPL (version 2) licence (see licence file for details).
//
// (c) R F L Evans 2014. All rights reserved.
//
//-----------------------------------------------------------------------------

// System headers
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

// Program headers
#include "atoms.hpp"
#include "errors.hpp"
#include "material.hpp"
#include "random.hpp"
#include "sim.hpp"
#include "sld.hpp"
#include "stats.hpp"
#include "vio.hpp"
#include "vmpi.hpp"
#include "program.hpp"

namespace{

   // marker for the integrator seed record appended to checkpoint files (ASCII "INTSEED1" in little endian)
   const uint64_t integration_seed_tag = 0x3144454553544E49;

   // integrator seed which started the loaded checkpoint chain, -1 if not recorded by an older checkpoint
   bool chain_seed_loaded = false;
   int32_t chain_seed = -1;

   // markers for the optional records appended after the seed record (ASCII "ALTFLD03" and "ENDCND01" in little endian)
   const uint64_t alternating_field_tag = 0x3330444C46544C41;
   const uint64_t end_condition_tag = 0x3130444E43444E45;

   // read the accumulated statistics, restoring them only if restore is true, otherwise the records are just consumed
   void load_statistics(std::ifstream& chkfile, const bool restore){
      stats::system_magnetization.load_checkpoint(chkfile,restore);
      stats::grain_magnetization.load_checkpoint(chkfile,restore);
      stats::material_magnetization.load_checkpoint(chkfile,restore);
      stats::material_grain_magnetization.load_checkpoint(chkfile,restore);
      stats::height_magnetization.load_checkpoint(chkfile,restore);
      stats::material_height_magnetization.load_checkpoint(chkfile,restore);
      stats::material_grain_height_magnetization.load_checkpoint(chkfile,restore);

      stats::system_specific_heat.load_checkpoint(chkfile,restore);
      stats::grain_specific_heat.load_checkpoint(chkfile,restore);
      stats::material_specific_heat.load_checkpoint(chkfile,restore);

      stats::system_susceptibility.load_checkpoint(chkfile,restore);
      stats::grain_susceptibility.load_checkpoint(chkfile,restore);
      stats::material_susceptibility.load_checkpoint(chkfile,restore);
   }

   // check that the input uses the same alternating field protocol and restore its schedule state
   void load_alternating_field_record(std::ifstream& chkfile){

      int32_t mode32, component32, method32;
      uint64_t interval_steps, buffer, partial_time;
      double trigger, dt;
      sim::alternating_field_state_t state;
      chkfile.read((char*)&mode32,sizeof(int32_t));
      chkfile.read((char*)&interval_steps,sizeof(uint64_t));
      chkfile.read((char*)&component32,sizeof(int32_t));
      chkfile.read((char*)&trigger,sizeof(double));
      chkfile.read((char*)&method32,sizeof(int32_t));
      chkfile.read((char*)&buffer,sizeof(uint64_t));
      chkfile.read((char*)&dt,sizeof(double));
      chkfile.read((char*)&partial_time,sizeof(uint64_t));
      chkfile.read((char*)&state.base_field,sizeof(double));
      chkfile.read((char*)&state.base_direction[0],sizeof(double)*3);
      chkfile.read((char*)&state.reversals,sizeof(uint64_t));
      chkfile.read((char*)&state.pending,sizeof(bool));
      chkfile.read((char*)&state.detection_time,sizeof(uint64_t));
      if(!chkfile) err::zexit("Incomplete alternating field record in checkpoint file");

      // the saved schedule state is only correct for the same protocol, so do not reinterpret it with changed input
      if(!sim::alternating_field_enabled) err::zexit("Checkpoint was saved with sim:alternating-applied-field, which must also be enabled to continue it");
      const bool same_protocol = mode32 == int32_t(sim::alternating_field_mode) && interval_steps == sim::alternating_field_interval_steps &&
                                 component32 == int32_t(sim::alternating_field_component) && trigger == sim::alternating_field_trigger &&
                                 method32 == int32_t(sim::alternating_field_trigger_method) && buffer == sim::alternating_field_buffer &&
                                 state.base_field == sim::H_applied && state.base_direction[0] == sim::H_vec[0] &&
                                 state.base_direction[1] == sim::H_vec[1] && state.base_direction[2] == sim::H_vec[2];
      if(!same_protocol) err::zexit("Alternating field settings or applied field differ from those saved in the checkpoint");

      // the schedule is in time steps, so a different time step changes its duration, and the trigger is only checked at output samples
      if(dt != mp::dt_SI) err::zexit("sim:time-step differs from the checkpoint, which would change the duration of the alternating field schedule");
      if(sim::alternating_field_mode == sim::alternating_field_trigger_mode && partial_time != sim::partial_time)
         err::zexit("sim:time-steps-increment differs from the checkpoint, which would change when the alternating field trigger is detected");

      state.switch_time = state.detection_time + sim::alternating_field_buffer;
      sim::alternating_field_state = state;

   }

   // read the latched end condition, which is only used when the input also enables an end condition
   void load_end_condition_record(std::ifstream& chkfile){
      bool triggered = false;
      uint64_t end_time = 0;
      chkfile.read((char*)&triggered,sizeof(bool));
      chkfile.read((char*)&end_time,sizeof(uint64_t));
      if(!chkfile) err::zexit("Incomplete end condition record in checkpoint file");
      sim::end_condition_triggered = triggered;
      sim::end_condition_end_time = end_time;
   }

}

//-----------------------------------------------------------------------------
// Function to save checkpoint file
//-----------------------------------------------------------------------------
void save_checkpoint(){

   // convert number of atoms, rank and time to standard long int
   uint64_t natoms64 = uint64_t(atoms::num_atoms-vmpi::num_halo_atoms);
   int64_t time64 = int64_t(sim::time);
   int64_t eqtime64 = int64_t(sim::equilibration_time);
   int64_t parity64 = int64_t(sim::parity);
   int64_t iH64 = int64_t(sim::iH);
   double temp = sim::temperature;
   int64_t output_atoms_file_counter64 = int64_t(sim::output_atoms_file_counter);
   int64_t output_cells_file_counter64 = int64_t(sim::output_cells_file_counter);
   int64_t output_rate_counter64 = int64_t(sim::output_rate_counter);
   double constr_theta = sim::constraint_theta;
   double constr_phi   = sim::constraint_phi;
   bool flag_constraint_theta_changed = sim::constraint_theta_changed;
   bool flag_constraint_phi_changed   = sim::constraint_phi_changed;

   // determine checkpoint file name
   std::stringstream chkfilenamess;
   chkfilenamess << "vampire" << vmpi::my_rank << ".chk";
   std::string chkfilename = chkfilenamess.str();

   // open checkpoint file
   std::ofstream chkfile;
   chkfile.open(chkfilename.c_str(),std::ios::binary);

   // check for open file
   if(!chkfile.is_open()){
      terminaltextcolor(RED);
      std::cerr << "Error: Unable to open checkpoint file " << chkfilename << " for writing. Exiting." << std::endl;
      terminaltextcolor(WHITE);
      zlog << zTs() << "Error: Unable to open checkpoint file " << chkfilename << " for writing. Exiting." << std::endl;
      err::vexit();
   }

   // get state of random number generator
   std::vector<uint32_t> mt_state(624); // 624 is hard coded in mt implementation. uint64 assumes same size as unsigned long
   int32_t mt_p=0; // position in rng state
   mt_p=mtrandom::grnd.get_state(mt_state);

   // write checkpoint variables to file
   chkfile.write(reinterpret_cast<const char*>(&natoms64),sizeof(uint64_t));
   chkfile.write(reinterpret_cast<const char*>(&time64),sizeof(int64_t));
   chkfile.write(reinterpret_cast<const char*>(&eqtime64),sizeof(int64_t));
   chkfile.write(reinterpret_cast<const char*>(&parity64),sizeof(int64_t));
   chkfile.write(reinterpret_cast<const char*>(&iH64),sizeof(int64_t));
   chkfile.write(reinterpret_cast<const char*>(&temp),sizeof(double));
   chkfile.write(reinterpret_cast<const char*>(&constr_theta),sizeof(double));
   chkfile.write(reinterpret_cast<const char*>(&constr_phi),sizeof(double));
   chkfile.write(reinterpret_cast<const char*>(&flag_constraint_theta_changed),sizeof(bool));
   chkfile.write(reinterpret_cast<const char*>(&flag_constraint_phi_changed  ),sizeof(bool));
   chkfile.write(reinterpret_cast<const char*>(&output_atoms_file_counter64),sizeof(int64_t));
   chkfile.write(reinterpret_cast<const char*>(&output_cells_file_counter64),sizeof(int64_t));
   chkfile.write(reinterpret_cast<const char*>(&output_rate_counter64),sizeof(int64_t));
   chkfile.write(reinterpret_cast<const char*>(&mt_p),sizeof(int32_t));
   chkfile.write(reinterpret_cast<const char*>(&mt_state[0]),sizeof(uint32_t)*mt_state.size());

   // write spin array to file
   chkfile.write(reinterpret_cast<const char*>(&atoms::x_spin_array[0]),sizeof(double)*natoms64);
   chkfile.write(reinterpret_cast<const char*>(&atoms::y_spin_array[0]),sizeof(double)*natoms64);
   chkfile.write(reinterpret_cast<const char*>(&atoms::z_spin_array[0]),sizeof(double)*natoms64);

   // Save positions and velocities only for SLD simulations
   if(sld::enabled){

      // write position array to file
      chkfile.write(reinterpret_cast<const char*>(&atoms::x_coord_array[0]),sizeof(double)*natoms64);
      chkfile.write(reinterpret_cast<const char*>(&atoms::y_coord_array[0]),sizeof(double)*natoms64);
      chkfile.write(reinterpret_cast<const char*>(&atoms::z_coord_array[0]),sizeof(double)*natoms64);

     // write velocity array to file
     chkfile.write(reinterpret_cast<const char*>(&atoms::x_velo_array[0]),sizeof(double)*natoms64);
     chkfile.write(reinterpret_cast<const char*>(&atoms::y_velo_array[0]),sizeof(double)*natoms64);
     chkfile.write(reinterpret_cast<const char*>(&atoms::z_velo_array[0]),sizeof(double)*natoms64);

  }

   // write statistical properties to file
   stats::system_magnetization.save_checkpoint(chkfile);
   stats::grain_magnetization.save_checkpoint(chkfile);
   stats::material_magnetization.save_checkpoint(chkfile);
   stats::material_grain_magnetization.save_checkpoint(chkfile);
   stats::height_magnetization.save_checkpoint(chkfile);
   stats::material_height_magnetization.save_checkpoint(chkfile);
   stats::material_grain_height_magnetization.save_checkpoint(chkfile);

   stats::system_specific_heat.save_checkpoint(chkfile);
   stats::grain_specific_heat.save_checkpoint(chkfile);
   stats::material_specific_heat.save_checkpoint(chkfile);

   stats::system_susceptibility.save_checkpoint(chkfile);
   stats::grain_susceptibility.save_checkpoint(chkfile);
   stats::material_susceptibility.save_checkpoint(chkfile);

   // append the seed which started this simulation after all other data so older versions can still read the file
   const int32_t seed32 = chain_seed_loaded ? chain_seed : int32_t(mtrandom::integration_seed);
   chkfile.write(reinterpret_cast<const char*>(&integration_seed_tag),sizeof(uint64_t));
   chkfile.write(reinterpret_cast<const char*>(&seed32),sizeof(int32_t));

   // append the alternating field protocol and schedule state, only when enabled so other checkpoints are unchanged
   if(sim::alternating_field_enabled){
      const int32_t mode32 = int32_t(sim::alternating_field_mode);
      const int32_t component32 = int32_t(sim::alternating_field_component);
      const int32_t method32 = int32_t(sim::alternating_field_trigger_method);
      const uint64_t partial_time = sim::partial_time;
      const sim::alternating_field_state_t& state = sim::alternating_field_state;
      chkfile.write(reinterpret_cast<const char*>(&alternating_field_tag),sizeof(uint64_t));
      chkfile.write(reinterpret_cast<const char*>(&mode32),sizeof(int32_t));
      chkfile.write(reinterpret_cast<const char*>(&sim::alternating_field_interval_steps),sizeof(uint64_t));
      chkfile.write(reinterpret_cast<const char*>(&component32),sizeof(int32_t));
      chkfile.write(reinterpret_cast<const char*>(&sim::alternating_field_trigger),sizeof(double));
      chkfile.write(reinterpret_cast<const char*>(&method32),sizeof(int32_t));
      chkfile.write(reinterpret_cast<const char*>(&sim::alternating_field_buffer),sizeof(uint64_t));
      chkfile.write(reinterpret_cast<const char*>(&mp::dt_SI),sizeof(double));
      chkfile.write(reinterpret_cast<const char*>(&partial_time),sizeof(uint64_t));
      chkfile.write(reinterpret_cast<const char*>(&state.base_field),sizeof(double));
      chkfile.write(reinterpret_cast<const char*>(&state.base_direction[0]),sizeof(double)*3);
      chkfile.write(reinterpret_cast<const char*>(&state.reversals),sizeof(uint64_t));
      chkfile.write(reinterpret_cast<const char*>(&state.pending),sizeof(bool));
      chkfile.write(reinterpret_cast<const char*>(&state.detection_time),sizeof(uint64_t));
   }

   // append the end condition latch, only when an end condition is enabled
   if(sim::end_condition_enabled){
      chkfile.write(reinterpret_cast<const char*>(&end_condition_tag),sizeof(uint64_t));
      chkfile.write(reinterpret_cast<const char*>(&sim::end_condition_triggered),sizeof(bool));
      chkfile.write(reinterpret_cast<const char*>(&sim::end_condition_end_time),sizeof(uint64_t));
   }

   sld::save_phonon_pump_checkpoint(chkfile);

   // close checkpoint file
   chkfile.close();

   // log writing checkpoint file (only for non-continuous checkpoint files)
   if(!sim::save_checkpoint_continuous_flag) zlog << zTs() << "Checkpoint file written to disk." << std::endl;

   return;

}

//-----------------------------------------------------------------------------
// Function to save checkpoint file
//-----------------------------------------------------------------------------
void load_checkpoint(){

   // convert number of atoms, rank and time to standard long int
   uint64_t natoms64;
   int64_t time64;
   int64_t eqtime64;
   int64_t parity64;
   int64_t iH64;
   double temp;
   int64_t output_atoms_file_counter64;
   int64_t output_cells_file_counter64;
   int64_t output_rate_counter64;
   double constr_theta;
   double constr_phi;
   bool flag_constraint_theta_changed;
   bool flag_constraint_phi_changed  ;

   // variables for loading state of random number generator
   std::vector<uint32_t> mt_state(624); // 624 is hard coded in mt implementation. uint64 assumes same size as unsigned long
   int32_t mt_p=0; // position in rng state

   // determine checkpoint file name
   std::stringstream chkfilenamess;
   chkfilenamess << "vampire" << vmpi::my_rank << ".chk";
   std::string chkfilename = chkfilenamess.str();

   // open checkpoint file
   std::ifstream chkfile;
   chkfile.open(chkfilename.c_str(),std::ios::binary);

   // check for open file
   if(!chkfile.is_open()){
      terminaltextcolor(RED);
      std::cerr << "Error: Unable to open checkpoint file " << chkfilename << " for reading. Exiting." << std::endl;
      std::cerr << "Info: sim:continue may be specified in the input file which requires a valid checkpoint file." << std::endl;
      terminaltextcolor(WHITE);
      zlog << zTs() << "Error: Unable to open checkpoint file " << chkfilename << " for reading. Exiting." << std::endl;
      zlog << zTs() << "Info: sim:continue may be specified in the input file which requires a valid checkpoint file." << std::endl;
      err::vexit();
   }

   // Set flag to true do determine that this is the beginning of the simulation
   sim::checkpoint_loaded_flag=true;
   zlog << zTs() << "Flag:checkpoint_loaded_flag = " << sim::checkpoint_loaded_flag <<std::endl;

   // read checkpoint variables from file
   chkfile.read((char*)&natoms64,sizeof(uint64_t));
   chkfile.read((char*)&time64,sizeof(int64_t));
   chkfile.read((char*)&eqtime64,sizeof(int64_t));
   chkfile.read((char*)&parity64,sizeof(int64_t));
   chkfile.read((char*)&iH64,sizeof(int64_t));
   chkfile.read((char*)&temp,sizeof(double));
   chkfile.read((char*)&constr_theta,sizeof(double));
   chkfile.read((char*)&constr_phi,sizeof(double));
   chkfile.read((char*)&flag_constraint_theta_changed,sizeof(bool));
   chkfile.read((char*)&flag_constraint_phi_changed  ,sizeof(bool));
   chkfile.read((char*)&output_atoms_file_counter64,sizeof(int64_t));
   chkfile.read((char*)&output_cells_file_counter64,sizeof(int64_t));
   chkfile.read((char*)&output_rate_counter64,sizeof(int64_t));
   chkfile.read((char*)&mt_p,sizeof(int32_t));
   chkfile.read((char*)&mt_state[0],sizeof(uint32_t)*mt_state.size());

   //std::cout << "random generator state loaded = " << mt_p << std::endl;
   // if continuing set state of rng
   if(sim::load_checkpoint_continue_flag) mtrandom::grnd.set_state(mt_state, mt_p);

   // check for rational number of atoms
   if(static_cast<uint64_t>(atoms::num_atoms-vmpi::num_halo_atoms) != natoms64){
      terminaltextcolor(RED);
      std::cerr << "Error: Mismatch between number of atoms in checkpoint file (" << natoms64 << ") and number of generated atoms (" << atoms::num_atoms-vmpi::num_halo_atoms << "). Exiting." << std::endl;
      terminaltextcolor(WHITE);
      zlog << zTs() << "Error: Mismatch between number of atoms in checkpoint file (" << natoms64 << ") and number of generated atoms (" << atoms::num_atoms-vmpi::num_halo_atoms << "). Exiting." << std::endl;
      err::vexit();
   }

   // Load saved parameters if simulation continuing
   if(sim::load_checkpoint_continue_flag){
      sim::parity = parity64;
      sim::iH = iH64;
      sim::time = time64;
      sim::equilibration_time = eqtime64;
      sim::temperature = temp;
      sim::output_atoms_file_counter = output_atoms_file_counter64;
      sim::output_cells_file_counter = output_cells_file_counter64;
      sim::output_rate_counter = output_rate_counter64;
      sim::constraint_theta = constr_theta;
      sim::constraint_phi = constr_phi;
      sim::constraint_theta_changed = flag_constraint_theta_changed;
      sim::constraint_phi_changed   = flag_constraint_phi_changed  ;
   }

   // Load spin positions
   chkfile.read((char*)&atoms::x_spin_array[0],sizeof(double)*natoms64);
   chkfile.read((char*)&atoms::y_spin_array[0],sizeof(double)*natoms64);
   chkfile.read((char*)&atoms::z_spin_array[0],sizeof(double)*natoms64);

   // Load positions and velocities only for SLD simulations
   if(sld::enabled){

      chkfile.read((char*)&atoms::x_coord_array[0],sizeof(double)*natoms64);
      chkfile.read((char*)&atoms::y_coord_array[0],sizeof(double)*natoms64);
      chkfile.read((char*)&atoms::z_coord_array[0],sizeof(double)*natoms64);

      chkfile.read((char*)&atoms::x_velo_array[0],sizeof(double)*natoms64);
      chkfile.read((char*)&atoms::y_velo_array[0],sizeof(double)*natoms64);
      chkfile.read((char*)&atoms::z_velo_array[0],sizeof(double)*natoms64);

   }

   // restore the SLED thermostat state
   if(sim::load_checkpoint_continue_flag) sld::restore_thermostat_checkpoint();

   // refresh restored spins and coordinates
   #ifdef MPICF
      if(sld::enabled){
         vmpi::mpi_init_halo_swap();
         vmpi::mpi_complete_halo_swap();
         vmpi::mpi_init_halo_swap_coords();
         vmpi::mpi_complete_halo_swap_coords();
      }
   #endif
   
   // The optional records after the statistics restore the alternating field state, which selects the field used by the first
   // instantaneous statistics. So consume the accumulated statistics without restoring them, read the optional records, calculate
   // the instantaneous statistics once, then return and restore the accumulated statistics, so that the checkpoint sample is not included twice.
   const std::streampos statistics_start = chkfile.tellg();
   load_statistics(chkfile, false);

   // read the optional seed record, which older checkpoint files end before
   uint64_t seed_tag = 0;
   int32_t seed32 = -1;
   const bool seed_record = chkfile.read((char*)&seed_tag,sizeof(uint64_t)) && seed_tag == integration_seed_tag && chkfile.read((char*)&seed32,sizeof(int32_t));
   if(!seed_record) seed32 = -1;

   // read the optional records which follow the seed record when continuing, a restart uses the settings in the input
   bool alternating_field_record_loaded = false;
   bool phonon_pump_record_loaded = false;
   if(seed_record && sim::load_checkpoint_continue_flag){
      uint64_t tag = 0;
      while(chkfile.read((char*)&tag,sizeof(uint64_t))){
         if(tag == alternating_field_tag){
            load_alternating_field_record(chkfile);
            alternating_field_record_loaded = true;
         }
         else if(tag == 0x3230504D55504E50){ // PNPUMP02
            sld::load_phonon_pump_checkpoint(chkfile);
            phonon_pump_record_loaded = true;
         }
         else if(tag == end_condition_tag) load_end_condition_record(chkfile);
         else break;
      }
   }
   if(sim::load_checkpoint_continue_flag && sim::alternating_field_enabled && !alternating_field_record_loaded)
      err::zexit("Cannot continue an alternating field from a checkpoint without alternating field state");

   if(sim::load_checkpoint_continue_flag) sld::check_phonon_pump_checkpoint(phonon_pump_record_loaded);

   // select the time series field for the restored stage and reversal count, then calculate the instantaneous statistics from the restored atom state
   if(sim::load_checkpoint_continue_flag && program::program == 1) program::set_time_series_field(true);
   stats::update();

   // restore the accumulated statistics
   chkfile.clear();
   chkfile.seekg(statistics_start);
   load_statistics(chkfile, sim::load_checkpoint_continue_flag);

   // close checkpoint file
   chkfile.close();

   // log reading checkpoint file
   zlog << zTs() << "Checkpoint file loaded at sim::time " << sim::time << "." << std::endl;

   // a continued simulation keeps the original seed for later checkpoints, a restart is seeded from the input file
   std::stringstream seedss;
   if(seed32 >= 0) seedss << seed32;
   else seedss << "not recorded (older checkpoint format)";
   if(sim::load_checkpoint_continue_flag){
      chain_seed_loaded = true;
      chain_seed = seed32;
      if(vmpi::master) std::cout << "Integrator random seed (from checkpoint): " << seedss.str() << std::endl;
      zlog << zTs() << "Integrator random seed (from checkpoint): " << seedss.str() << std::endl;
   }
   else{
      zlog << zTs() << "Integrator random seed of loaded checkpoint: " << seedss.str() << ", restarting with integrator random seed " << mtrandom::integration_seed << std::endl;
   }

   return;

}
