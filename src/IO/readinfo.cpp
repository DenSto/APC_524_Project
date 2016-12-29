#include <stdio.h>
#include <assert.h>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <string.h>

#include "input.hpp"
#include "../globals.hpp"
#include "../../libconfig-1.5/lib/libconfig.h++"

using namespace std;
using namespace libconfig;

int Input::readinfo(char *fname){

    Config cfg;
    // Read the file. If there is an error, report it and exit.
    try
    {
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

    try
    {
      const Setting &nCell = cfg.lookup("domain.nCell");
      if(nCell.getLength() == 3) {
        input_info_->nCell[0] = nCell[0];
        input_info_->nCell[1] = nCell[1];
        input_info_->nCell[2] = nCell[2];
      } else {
        input_info_->nCell[0] = nCell[0];
        input_info_->nCell[1] = nCell[0];
        input_info_->nCell[2] = nCell[0];
        cerr << "Error: nCell is not a 3 element array in input file."
           << endl << "Assuming nCell[0]=nCell[1]=nCell[2]=" 
           << input_info_->nCell[0] << "." << endl;
      }
    }
    catch(const SettingNotFoundException &nfex)
    {
      cerr << "Error: nCell not set in input file" << endl; 
      return(EXIT_FAILURE);
    }

// dummy code to be replaced ///////////////////////////////
    input_info_->xyz0[0]=  0.0;
    input_info_->xyz0[1]=  0.0;
    input_info_->xyz0[2]=  0.0;

    input_info_->Lxyz[0] = 1.5;
    input_info_->Lxyz[1] = 2.0;
    input_info_->Lxyz[2] = 2.5;
///////////////////////////////////////////////////////////

#if USE_MPI
    try
    {
      const Setting &nProc = cfg.lookup("domain.nProc");
      if(nProc.getLength() == 3) {
        input_info_->nProc[0] = nProc[0];
        input_info_->nProc[1] = nProc[1];
        input_info_->nProc[2] = nProc[2];
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
      cerr << "startTime not set in input file... " 
           << "Assuming startTime = 0" << endl;
      input_info_->t0 = 0.;
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
      input_info_->restart = cfg.lookup("initialization.restart");
    }
    catch(const SettingNotFoundException &nfex)
    {
      cerr << "restart not set in input file... " 
           << "Assuming restart = 0" << endl;
      input_info_->restart = 0.;
    }

    int nspecies;
    try
    {
      nspecies = cfg.lookup("initialization.particles.nSpecies");
      input_info_->nspecies = nspecies;
      if(nspecies < 1) {
        cerr << "Error: no species in simulation ..." << endl;
        return(EXIT_FAILURE); 
      }
    }
    catch(const SettingNotFoundException &nfex)
    {
      cerr << "Error: nSpecies not set in input file" << endl; 
      return(EXIT_FAILURE);
    }

    // allocate input_info
    mallocinfo(nspecies);
    //cerr << "mallocinfo in readinfo.cpp" << endl;

    try
    {
      const Setting &mass_ratio = 
           cfg.lookup("initialization.particles.mass_ratio");
      for(int i=0;i<nspecies;i++){
        input_info_->mass_ratio[i] = mass_ratio[i];
      }
      if(mass_ratio.getLength() > nspecies) {
        cerr << "Error: More mass ratios are specified than nspecies!" << endl;
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
      for(int i=0;i<nspecies;i++){
        input_info_->charge_ratio[i] = charge_ratio[i];
      }
      if(charge_ratio.getLength() > nspecies) {
        cerr << "Error: More charge ratios are specified than nspecies!" << endl;
      }
    }
    catch(const SettingNotFoundException &nfex)
    {
      cerr << "Error: charge_ratio not set in input file" << endl; 
      return(EXIT_FAILURE);
    }

    try
    {
      const Setting &dens_frac = 
           cfg.lookup("initialization.particles.dens_frac");
      for(int i=0;i<nspecies;i++){
        input_info_->dens_frac[i] = dens_frac[i];
      }
      if(dens_frac.getLength() > nspecies) {
        cerr << "Error: More fractional densities are specified than nspecies!" << endl;
      }
    }
    catch(const SettingNotFoundException &nfex)
    {
      cerr << "Error: charge_ratio not set in input file" << endl; 
      return(EXIT_FAILURE);
    }

    try
    {
      const Setting &nParticles = 
           cfg.lookup("initialization.particles.nParticles");
      if(nParticles.getType() == Setting::TypeInt) {
        int np = nParticles;
        input_info_->np = (long) np;
        if(input_info_->np < 0) {
          cerr << "Error: integer overflow..." 
               << "Use nParticles = #######L in input file"
               << "for long format or use scientific notation." << endl;
          return(EXIT_FAILURE); 
        }
      } else if (nParticles.getType() == Setting::TypeInt64) {
        input_info_->np = nParticles;
      } else if (nParticles.getType() == Setting::TypeFloat) {
        double np = nParticles;
        input_info_->np = (long) np;
      }
    }
    catch(const SettingNotFoundException &nfex)
    {
      cerr << "Error: nParticles not set in input file" << endl; 
      return(EXIT_FAILURE);
    }

// dummy code to be replaced //////////////////////////////////
    //input_info->dens    = 0.2;
    //input_info->temp    = 1.5;

    sprintf(input_info_->distname,"distribution.dat");

    // MPI can only Bcast C strings
    char (*parts_bound)[32] = input_info_->parts_bound;
    strcpy(parts_bound[0], "periodic"); // x -> Left
    strcpy(parts_bound[1], "periodic"); // x -> Right
    strcpy(parts_bound[2], "periodic"); // y -> Left
    strcpy(parts_bound[3], "periodic"); // y -> Right
    strcpy(parts_bound[4], "periodic"); // z -> Left
    strcpy(parts_bound[5], "periodic"); // z -> Right
 
    char (*fields_bound)[32] = input_info_->fields_bound;
    strcpy(fields_bound[0], "periodic"); // x -> Left
    strcpy(fields_bound[1], "periodic"); // x -> Right
    strcpy(fields_bound[2], "periodic"); // y -> Left
    strcpy(fields_bound[3], "periodic"); // y -> Right
    strcpy(fields_bound[4], "periodic"); // z -> Left
    strcpy(fields_bound[5], "periodic"); // z -> Right
////////////////////////////////////////////////////////////////

    return 0;
}

