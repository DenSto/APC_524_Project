#include <stdio.h>
#include <assert.h>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <string.h>

#include "input.hpp"
#include "../globals.hpp"
#include "../../libconfig-1.5/lib/libconfig.h++"

#define TRUE 1
#define FALSE 0

using namespace std;
using namespace libconfig;

int Input::readinfo(char *fname){

    Config cfg;
    // Read the file. If there is an error, report it and exit.

    cfg.setAutoConvert(true);

    try
    {
      cout << "Reading input file " << fname << endl;
      cfg.readFile(fname);
    }
    catch(const FileIOException &fioex)
    {
      cerr << "I/O error while reading file." << endl;
      return(EXIT_FAILURE);
    }
    catch(const ParseException &pex)
    {
      cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
                << " - " << pex.getError() << endl;
      return(EXIT_FAILURE);
    }

    /* Read domain inputs *********************************/
    try
    {
      const Setting &nCell = cfg.lookup("domain.nCell");
      if(nCell.getLength() == NDIM) {
        for(int i=0;i<NDIM;i++){
          input_info_->nCell[i] = nCell[i];
        }
      } else {
        for(int i=0;i<NDIM;i++){
          input_info_->nCell[i] = nCell[0];
        }
        cerr << "Caution: nCell is not a " << NDIM << " element array in input file!"
           << endl << "  Assuming nCell[i]=" 
           << input_info_->nCell[0] << "." << endl;
      }
    }
    catch(const SettingNotFoundException &nfex)
    {
      cerr << "Error: nCell not set in input file!" << endl; 
      return(EXIT_FAILURE);
    }

    try
    {
      const Setting &xyz0 = cfg.lookup("domain.xyz0");
      if(xyz0.getLength() == NDIM) {
        for(int i=0;i<NDIM;i++){
          input_info_->xyz0[i] = xyz0[i];
        }
      } else {
        for(int i=0;i<NDIM;i++){
          input_info_->xyz0[i] = xyz0[0];
        }
        cerr << "Caution: xyz0 is not a " << NDIM << " element array in input file."
           << endl << "Assuming xyz0[i]=" 
           << input_info_->xyz0[0] << "." << endl;
      }
    }
    catch(const SettingNotFoundException &nfex)
    {
      cerr << "Caution: xyz0 not set in input file!" << endl
           << "  Assuming xyz0[i]=0.0" << endl; 
    }

    try
    {
      const Setting &Lxyz = cfg.lookup("domain.Lxyz");
      if(Lxyz.getLength() == NDIM) {
        for(int i=0;i<NDIM;i++){
          input_info_->Lxyz[i] = Lxyz[i];
        }
      } else {
        for(int i=0;i<NDIM;i++){
          input_info_->Lxyz[i] = Lxyz[0];
        }
        cerr << "Caution: Lxyz is not a " << NDIM << " element array in input file."
           << endl << "Assuming Lxyz[i]=" 
           << input_info_->Lxyz[0] << "." << endl;
      }
    }
    catch(const SettingNotFoundException &nfex)
    {
      cerr << "Error: Lxyz not set in input file!" << endl;
      return(EXIT_FAILURE);
    }


#if USE_MPI
    try
    {
      const Setting &nProc = cfg.lookup("domain.nProc");
      if(nProc.getLength() == NDIM) {
        for(int i=0;i<NDIM;i++){
          input_info_->nProc[i] = nProc[i];
        }
      } else if(nProc.getLength()==1) {
        input_info_->nProc[0] = nProc[0];
        input_info_->nProc[1] = 1;
        input_info_->nProc[2] = 1;
        cerr << "Error: nProc is not a 3 element array in input file."
           << endl << "Assuming nProc[0]=" << input_info_->nProc[0] 
           << ", nProc[1]=nProc[2]=1." << endl;
      } else {
        cerr << "Error: unrecognized nProc input format. Use" << endl
             << "nProc = [# # #]." << endl;
        return(EXIT_FAILURE);
      }
      if(input_info_->nProc[0]*input_info_->nProc[1]*input_info_->nProc[2]!=size_MPI)
      {
        cerr << "Error: nProc layout specified in input file does not match "
             << "number of MPI processes requested." << endl;
        return(EXIT_FAILURE);
      }
    }
    catch(const SettingNotFoundException &nfex)
    {
      cerr << "Error: nProc not set in input file" << endl; 
      return(EXIT_FAILURE);
    }
#else
    // serial case
    input_info_->nProc[0] = input_info_->nProc[1] = input_info_->nProc[2] = 1;
#endif 

    /* Read runtime inputs ********************************/
    try
    {
      input_info_->nt = cfg.lookup("runtime.nTimesteps");
    }
    catch(const SettingNotFoundException &nfex)
    {
      cerr << "Error: nTimesteps not set in input file" << endl; 
      return(EXIT_FAILURE);
    }

    try
    {
      input_info_->t0 = cfg.lookup("runtime.startTime");
    }
    catch(const SettingNotFoundException &nfex)
    {
      cerr << "Caution: startTime not set in input file!" << endl 
           << "  Assuming startTime = 0" << endl;
      input_info_->t0 = 0.;
    }

    try{input_info_->dens_phys = cfg.lookup("initialization.particles.dens_phys");
    }catch(const SettingNotFoundException &nfex){
      cerr << "Caution: dens_phys not set in input file!" << endl 
           << "  Assuming dens_phys = actual particle number density" << endl;
      input_info_->dens_phys = 0.0;
    }

    try
    {
      input_info_->debug = cfg.lookup("runtime.debug");
    }
    catch(const SettingNotFoundException &nfex)
    {
      input_info_->debug = 0;
    }

    try
    {
      input_info_->nstep_sort = cfg.lookup("runtime.nstep_sort");
    }
    catch(const SettingNotFoundException &nfex)
    {
      cerr << "nstep_sort not set. Using a sort cadence of 100 steps." << endl;
      input_info_->nstep_sort = 100;
    }
    
    int restart = 0;
    try
    {
      restart = cfg.lookup("initialization.restart");
      input_info_->restart = restart;
    }
    catch(const SettingNotFoundException &nfex)
    {
      cerr << "restart not set in input file... " 
           << "Assuming restart = 0" << endl;
      input_info_->restart = 0;
    }


    /* Read particle inputs *******************************/
    try
    {
      input_info_->relativity = cfg.lookup("initialization.particles.relativity");
    }
    catch(const SettingNotFoundException *nfex)
    {
      cerr << "relativity not set in input file..."
           << "Assuming relativity = 0" << endl;
      input_info_->relativity = 0;
    }
      

    int nspecies;
    try
    {
      nspecies = cfg.lookup("initialization.particles.nspecies");
      input_info_->nspecies = nspecies;
      if(nspecies < 1) {
        cerr << "Error: no species in simulation ..." << endl;
        return(EXIT_FAILURE); 
      }
      if(nspecies > NSPEC) {
        cerr << "Error: More species are specified than NSPEC!" << endl
             << "Change NSPCE to reserve more memory and recompile!" << endl;
        return(EXIT_FAILURE); 
      }
    }
    catch(const SettingNotFoundException &nfex)
    {
      cerr << "Error: nSpecies not set in input file" << endl; 
      return(EXIT_FAILURE);
    }

    try
    {
      const Setting &mass_ratio = 
           cfg.lookup("initialization.particles.mass_ratio");
      if(mass_ratio.getLength() != nspecies && !(nspecies==1 && mass_ratio.getLength()==0)) {
        cerr << "Error: length of mass ratios does not match nspecies!" << endl;
        return(EXIT_FAILURE);
      } else {
        if(nspecies==1 && mass_ratio.getLength()==0) input_info_->mass_ratio[0] = mass_ratio;
        else {
          for(int i=0;i<nspecies;i++){
            input_info_->mass_ratio[i] = mass_ratio[i];
          }
        }
      }
    }
    catch(const SettingNotFoundException &nfex)
    {
      cerr << "Error: mass_ratio not set in input file" << endl; 
      return(EXIT_FAILURE);
    }

    try
    {
      const Setting &charge_ratio = 
           cfg.lookup("initialization.particles.charge_ratio");
      if(charge_ratio.getLength() != nspecies && !(nspecies==1 && charge_ratio.getLength()==0)) {
        cerr << "Error: length of charge ratios does not match nspecies!" << endl;
        return(EXIT_FAILURE);
      } else {
        if(nspecies==1 && charge_ratio.getLength()==0) input_info_->charge_ratio[0] = charge_ratio;
        else {
          for(int i=0;i<nspecies;i++){
            input_info_->charge_ratio[i] = charge_ratio[i];
          }
        }
      }
    }
    catch(const SettingNotFoundException &nfex)
    {
      cerr << "Error: charge_ratio not set in input file" << endl; 
      return(EXIT_FAILURE);
    }
    try
    {
      const Setting &temperature = 
           cfg.lookup("initialization.particles.temperature");
      if(temperature.getLength() != nspecies && !(nspecies==1 && temperature.getLength()==0)) {
        cerr << "Error: length of temp does not match nspecies!" << endl;
        return(EXIT_FAILURE);
      } else {
        if(nspecies==1 && temperature.getLength()==0) input_info_->temp[0] = temperature;
        else {
          for(int i=0;i<nspecies;i++){
            input_info_->temp[i] = temperature[i];
          }
        }
      }
    }
    catch(const SettingNotFoundException &nfex)
    {
      cerr << "Error: temperature not set in input file" << endl; 
      return(EXIT_FAILURE);
    }

    try
    {
      const Setting &dens_frac = 
           cfg.lookup("initialization.particles.dens_frac");
      if(dens_frac.getLength() != nspecies && !(nspecies==1 && dens_frac.getLength()==0)) {
        cerr << "Error: length of dens_frac does not match nspecies!" << endl;
        return(EXIT_FAILURE);
      } else {
        if(nspecies==1 && dens_frac.getLength()==0) input_info_->dens_frac[0] = dens_frac;
        else {
          for(int i=0;i<nspecies;i++){
            input_info_->dens_frac[i] = dens_frac[i];
          }
        }
      }
    }
    catch(const SettingNotFoundException &nfex)
    {
      if(nspecies==1) input_info_->dens_frac[0] = 1;
      else {
        cerr << "Error: dens_frac not set in input file" << endl; 
        return(EXIT_FAILURE);
      }
    }

    try
    {
      const Setting &isTestParticle = 
           cfg.lookup("initialization.particles.isTestParticle");
      if(isTestParticle.getLength() != nspecies && !(nspecies==1 && isTestParticle.getLength()==0)) {
        cerr << "Error: length of isTestParticle does not match nspecies!" << endl;
        return(EXIT_FAILURE);
      } else {
        if(nspecies==1 && isTestParticle.getLength()==0) input_info_->isTestParticle[0] = isTestParticle;
        else {
          for(int i=0;i<nspecies;i++){
            input_info_->isTestParticle[i] = isTestParticle[i];
          }
        }
      }
    }
    catch(const SettingNotFoundException &nfex)
    {
	cerr << "Warning: isTestParticle not set. Assuming no test particle species." << endl;
        for(int i=0;i<nspecies;i++){
          input_info_->isTestParticle[i] = 0;
        }
	
    }

    try
    {
      const Setting &nParticles_tot = 
           cfg.lookup("initialization.particles.nParticles_tot");
      if(nParticles_tot.getType() == Setting::TypeInt) {
        int nparts_tot = nParticles_tot;
        input_info_->nparticles_tot = (long) nparts_tot;
        if(input_info_->nparticles_tot < 0) {
          cerr << "Error: integer overflow..." 
               << "Use nParticles = #######L in input file"
               << "for long format or use scientific notation." << endl;
          return(EXIT_FAILURE); 
        }
      } else if (nParticles_tot.getType() == Setting::TypeInt64) {
        input_info_->nparticles_tot = nParticles_tot;
      } else if (nParticles_tot.getType() == Setting::TypeFloat) {
        double nparts_tot = nParticles_tot;
        input_info_->nparticles_tot = (long) nparts_tot;
      }
    }
    catch(const SettingNotFoundException &nfex)
    {
      cerr << "Error: nParticles_tot not set in input file" << endl; 
      return(EXIT_FAILURE);
    }
    
    char (*parts_bound)[NCHAR] = input_info_->parts_bound;
    try
    {
      const Setting &partbc = cfg.lookup("boundary.particles.conditions"); 
      //cerr << "Error: List is of type " << Setting::TypeList << endl; 
      //cerr << "Error: Group is of type " << Setting::TypeGroup << endl; 
      //cerr << "Error: Array is of type " << Setting::TypeArray << endl; 
      if(partbc.getType() == Setting::TypeArray) {
        if(partbc.getLength() != 2*NDIM) {
          cerr << "Error: length of particle boundary condition is" 
               << partbc.getLength() << "." << endl
               << "This does not match the number of dimensions!"  << endl;
          return(EXIT_FAILURE);
        } else {
          for(int i=0;i<2*NDIM;i++){
            if(partbc[i].getType() == Setting::TypeString){
              strcpy(parts_bound[i], partbc[i]); 
            } else {
              cerr << "Error: each element of particle boundary should be a string!" << endl; 
              return(EXIT_FAILURE);
            }
          }
        }
      } else {
          cerr << "Error: particle boundary conditions should be array of strings!" 
               << endl << "  right now it is of Type " << partbc.getType() << endl; 
          return(EXIT_FAILURE);
      }
    }
    catch(const SettingNotFoundException &nfex)
    {
      cerr << "Caution: particle boundary not set in input file!" << endl
           << "  Assume periodic boundary conditions." << endl; 
      for(int i=0;i<2*NDIM;i++){
        strcpy(parts_bound[i], "periodic"); 
      }
    }

    /* Read fields inputs *********************************/
    try{input_info_->electrostatic = cfg.lookup("initialization.fields.electrostatic");
    }catch(const SettingNotFoundException *nfex){
      cerr << "electrostatic not set in input file..."
           << "Assuming electrostatic = 0" << endl;
      input_info_->electrostatic = 0;
    }

    try {
      const Setting &fieldinit = cfg.lookup("initialization.fields.init"); 
      if(fieldinit.getType() == Setting::TypeString){
        strcpy(input_info_->fields_init, fieldinit); 
      } else {
        cerr << "Error: field initialization method should be a string!" << endl; 
        return(EXIT_FAILURE);
      }
    } catch(const SettingNotFoundException &nfex){
      cerr << "Error: field initialization method not set in input file!" << endl;
      return(EXIT_FAILURE);
    }

    // read constant fields
    try{const Setting &B0 = cfg.lookup("initialization.fields.B0");
      if(B0.getLength() == NDIM) {
        for(int i=0;i<NDIM;i++){input_info_->B0[i] = B0[i];}
      } else {
        for(int i=0;i<NDIM;i++){input_info_->B0[i] = 0.0;}
        cerr << "Caution: B0 is not a " << NDIM << " element array in input file."
           << endl << "Assuming B0[i] = 0.0" << endl;
      }
    }catch(const SettingNotFoundException &nfex){
      for(int i=0;i<NDIM;i++){input_info_->B0[i] = 0.0;}
      cerr << "Caution: B0 not set in input file!" << endl
           << "  Assuming B0[i]=0.0" << endl; 
    }
   
    try{const Setting &E0 = cfg.lookup("initialization.fields.E0");
      if(E0.getLength() == NDIM) {
        for(int i=0;i<NDIM;i++){input_info_->E0[i] = E0[i];}
      } else {
        for(int i=0;i<NDIM;i++){input_info_->E0[i] = 0.0;}
        cerr << "Caution: E0 is not a " << NDIM << " element array in input file."
           << endl << "Assuming E0[i] = 0.0" << endl;
      }
    }catch(const SettingNotFoundException &nfex){
      for(int i=0;i<NDIM;i++){input_info_->B0[i] = 0.0;}
      cerr << "Caution: E0 not set in input file!" << endl
           << "  Assuming E0[i]=0.0" << endl; 
    }
    

    char (*fields_bound)[NCHAR] = input_info_->fields_bound;
    int isexternal = 0;
    try
    {
      const Setting &fieldbc = cfg.lookup("boundary.fields.conditions"); 
      //cerr << "Error: List is of type " << Setting::TypeList << endl; 
      //cerr << "Error: Group is of type " << Setting::TypeGroup << endl; 
      //cerr << "Error: Array is of type " << Setting::TypeArray << endl; 
      if(fieldbc.getType() == Setting::TypeArray) {
        if(fieldbc.getLength() != 2*NDIM) {
          cerr << "Error: length of field boundary condition is" 
               << fieldbc.getLength() << "." << endl
               << "This does not match the number of dimensions!"  << endl;
          return(EXIT_FAILURE);
        } else {
          for(int i=0;i<2*NDIM;i++){
            if(fieldbc[i].getType() == Setting::TypeString){
              strcpy(fields_bound[i], fieldbc[i]); 
              if(strcmp(fieldbc[i],"external")==0){isexternal+=1;}
            } else {
              cerr << "Error: each element of field boundary should be a string!" << endl; 
              return(EXIT_FAILURE);
            }
          }
        }
      } else {
          cerr << "Error: field boundary conditions should be array of strings!" 
               << endl << "  right now it is of Type " << fieldbc.getType() << endl; 
          return(EXIT_FAILURE);
      }
    }
    catch(const SettingNotFoundException &nfex)
    {
      cerr << "Caution: field boundary not set in input file!" << endl
           << "  Assume periodic boundary conditions." << endl; 
      for(int i=0;i<2*NDIM;i++){
        strcpy(fields_bound[i], "periodic"); 
      }
    }
    
    if(isexternal>0){
      cout << "    There are " << isexternal << " external field boundaries" << endl;
      // reading external fields parameters    

      int nwaves = 0;
      try
      {
        nwaves = cfg.lookup("boundary.fields.external.nwaves");
        input_info_->nwaves = nwaves;
        if(nwaves < 1) {
          cerr << "Error: no external field specified ..." << endl;
          return(EXIT_FAILURE); 
        }
        if(nwaves > NWAVE) {
          cerr << "Error: More external waves are specified than NWAVE!" << endl
               << "Change NWAVE to reserve more memory and recompile!" << endl;
          return(EXIT_FAILURE); 
        }
      }
      catch(const SettingNotFoundException &nfex)
      {
        cerr << "Caution: nwave not set in input file" << endl
             << "  Assume zero boundary conditions." << endl; 
        input_info_->nwaves = 0;
      }
    
      try
      {
        const Setting &inSide = 
             cfg.lookup("boundary.fields.external.inSide");
        if(inSide.getLength() != nwaves) {
          cerr << "Error: length of inSide does not match nwaves!" << endl;
          return(EXIT_FAILURE);
        } else {
          for(int i=0;i<nwaves;i++){
            input_info_->inSide[i] = inSide[i];
          }
        }
      }
      catch(const SettingNotFoundException &nfex)
      {
        cerr << "Error: inSide not set in input file" << endl; 
        return(EXIT_FAILURE);
      }

      try
      {
        const Setting &inPolE = 
             cfg.lookup("boundary.fields.external.inPolE");
        if(inPolE.getLength() != nwaves) {
          cerr << "Error: length of inPol does not match nwaves!" << endl;
          return(EXIT_FAILURE);
        } else {
          for(int i=0;i<nwaves;i++){
            input_info_->inPolE[i] = inPolE[i];
          }
        }
      }
      catch(const SettingNotFoundException &nfex)
      {
        cerr << "Error: inPolE not set in input file" << endl; 
        return(EXIT_FAILURE);
      }

      try{const Setting &peakamps= cfg.lookup("boundary.fields.external.peakamps");
        if(peakamps.getLength() != nwaves) {
          cerr << "Error: length of peakamps does not match nwaves!" << endl;
          return(EXIT_FAILURE);
        } else {
          for(int i=0;i<nwaves;i++){input_info_->peakamps[i] = peakamps[i];}
        }
      }catch(const SettingNotFoundException &nfex){
        cerr << "Error: peakamps not set in input file" << endl; 
        return(EXIT_FAILURE);
      }

      try{const Setting &omegas= cfg.lookup("boundary.fields.external.omegas");
        if(omegas.getLength() != nwaves) {
          cerr << "Error: length of omegas does not match nwaves!" << endl;
          return(EXIT_FAILURE);
        } else {
          for(int i=0;i<nwaves;i++){input_info_->omegas[i] = omegas[i];}
        }
      }catch(const SettingNotFoundException &nfex){
        cerr << "Error: omegas not set in input file" << endl; 
        return(EXIT_FAILURE);
      }

      try{const Setting &phases= cfg.lookup("boundary.fields.external.phases");
        if(phases.getLength() != nwaves) {
          cerr << "Error: length of phases does not match nwaves!" << endl;
          return(EXIT_FAILURE);
        } else {
          for(int i=0;i<nwaves;i++){input_info_->phases[i] = phases[i];}
        }
      }catch(const SettingNotFoundException &nfex){
        cerr << "Error: phases not set in input file" << endl; 
        return(EXIT_FAILURE);
      }

      try{const Setting &invWidths = cfg.lookup("boundary.fields.external.invWidths");
        if(invWidths.getLength() != nwaves) {
          cerr << "Error: length of invWidths does not match nwaves!" << endl;
          return(EXIT_FAILURE);
        } else {
          for(int i=0;i<nwaves;i++){input_info_->invWidths[i] = invWidths[i];}
        }
      }catch(const SettingNotFoundException &nfex){
        cerr << "Error: invWidths not set in input file" << endl; 
        return(EXIT_FAILURE);
      }
      
      try{const Setting &delays = cfg.lookup("boundary.fields.external.delays");
        if(delays.getLength() != nwaves) {
          cerr << "Error: length of delays does not match nwaves!" << endl;
          return(EXIT_FAILURE);
        } else {
          for(int i=0;i<nwaves;i++){input_info_->delays[i] = delays[i];}
        }
      }catch(const SettingNotFoundException &nfex){
        cerr << "Error: delays not set in input file" << endl; 
        return(EXIT_FAILURE);
      }

      // reading boundary conditions for Poisson solver to be used at initialiation step
      if(restart==0&&strcmp(input_info_->fields_init,"poisson")==0){
         try{const Setting &bound_phi=cfg.lookup("boundary.fields.external.poisson.bound_phi");
           if(bound_phi.getLength() != 6) {
             cerr << "Error: length of bound_phi is not 6!" << endl;
             return(EXIT_FAILURE);
           } else {
             for(int i=0;i<6;i++){input_info_->bound_phi[i] = bound_phi[i];}
           }
         }catch(const SettingNotFoundException &nfex){
           cerr << "Error: bound_phi not set in input file" << endl; 
           return(EXIT_FAILURE);
         }

         try{const Setting &bound_Ax=cfg.lookup("boundary.fields.external.poisson.bound_Ax");
           if(bound_Ax.getLength() != 6) {
             cerr << "Error: length of bound_Ax is not 6!" << endl;
             return(EXIT_FAILURE);
           } else {
             for(int i=0;i<6;i++){input_info_->bound_Ax[i] = bound_Ax[i];}
           }
         }catch(const SettingNotFoundException &nfex){
           cerr << "Error: bound_Ax not set in input file" << endl; 
           return(EXIT_FAILURE);
         }

         try{const Setting &bound_Ay=cfg.lookup("boundary.fields.external.poisson.bound_Ay");
           if(bound_Ay.getLength() != 6) {
             cerr << "Error: length of bound_Ay is not 6!" << endl;
             return(EXIT_FAILURE);
           } else {
             for(int i=0;i<6;i++){input_info_->bound_Ay[i] = bound_Ay[i];}
           }
         }catch(const SettingNotFoundException &nfex){
           cerr << "Error: bound_Ay not set in input file" << endl; 
           return(EXIT_FAILURE);
         }

         try{const Setting &bound_Az=cfg.lookup("boundary.fields.external.poisson.bound_Az");
           if(bound_Az.getLength() != 6) {
             cerr << "Error: length of bound_Az is not 6!" << endl;
             return(EXIT_FAILURE);
           } else {
             for(int i=0;i<6;i++){input_info_->bound_Az[i] = bound_Az[i];}
           }
         }catch(const SettingNotFoundException &nfex){
           cerr << "Error: bound_Az not set in input file" << endl; 
           return(EXIT_FAILURE);
         }
      }
    }else{
      input_info_->nwaves = 0;
    }
  
    /* Diagnostic inputs **********************************/
    try{input_info_->nstep_fields = cfg.lookup("diagnostics.nstep_fields");
    }catch(const SettingNotFoundException &nfex){
      cerr << "nstep_fields not set in input file..."
           << "Assuming nstep_fields = 10" << endl;
      input_info_->nstep_fields = 10;
    }

    try{input_info_->nstep_parts = cfg.lookup("diagnostics.nstep_parts");
    }catch(const SettingNotFoundException &nfex){
      cerr << "nstep_parts not set in input file..."
           << "Assuming nstep_parts = 10" << endl;
      input_info_->nstep_parts = 10;
    }

    try{input_info_->nstep_restart = cfg.lookup("diagnostics.nstep_restart");
    }catch(const SettingNotFoundException &nfex){
      cerr << "nwrite not set in input file..."
           << "Assuming nstep_restart = 100" << endl;
      input_info_->nstep_restart = 100;
    }

    // particle
    try{input_info_->output_pCount = cfg.lookup("diagnostics.particles.output_pCount");
    }catch(const SettingNotFoundException *nfex){
      cerr << "output_pCound not set in input file..."
           << "Assuming output_pCount = 0" << endl;
      input_info_->output_pCount = 0; // no tracking
    }
  
    //fields
    try {input_info_->which_fields = cfg.lookup("diagnostics.fields.which_fields");
    }catch(const SettingNotFoundException &nfex){
      cerr << "which_fields not set in input file..."
           << "Assuming which_fields = -1" << endl;
      input_info_->which_fields = -1; // no output
    }


////////////////////////////////////////////////////////////////

    return 0;
}

