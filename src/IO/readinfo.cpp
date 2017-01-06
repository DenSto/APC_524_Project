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
      if(nCell.getLength() == NDIM) {
        for(int i=0;i<NDIM;i++){
          input_info_->nCell[i] = nCell[i];
        }
      } else {
        for(int i=0;i<nCell.getLength();i++){
          input_info_->nCell[i] = nCell[0];
        }
        cerr << "Error: nCell is not a " << NDIM << " element array in input file."
           << endl << "Assuming nCell[i]=" 
           << input_info_->nCell[0] << "." << endl;
      }
    }
    catch(const SettingNotFoundException &nfex)
    {
      cerr << "Error: nCell not set in input file" << endl; 
      return(EXIT_FAILURE);
    }

// dummy code to be replaced ///////////////////////////////
    input_info_->relativity = 0;

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
      input_info_->restart = 0;
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
      if(mass_ratio.getLength() != nspecies) {
        cerr << "Error: length of mass ratios does not match nspecies!" << endl;
        return(EXIT_FAILURE);
      } else {
        for(int i=0;i<nspecies;i++){
          input_info_->mass_ratio[i] = mass_ratio[i];
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
      if(charge_ratio.getLength() != nspecies) {
        cerr << "Error: length of charge ratios does not match nspecies!" << endl;
        return(EXIT_FAILURE);
      } else {
        for(int i=0;i<nspecies;i++){
          input_info_->charge_ratio[i] = charge_ratio[i];
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
      const Setting &dens_frac = 
           cfg.lookup("initialization.particles.dens_frac");
      if(dens_frac.getLength() != nspecies) {
        cerr << "Error: length of fractional densities does not match nspecies!" << endl;
        return(EXIT_FAILURE);
      } else {
        for(int i=0;i<nspecies;i++){
          input_info_->dens_frac[i] = dens_frac[i];
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

    // diagnostics parameters
    try
    {
      input_info_->nwrite = cfg.lookup("diagnostics.nwrite");
    }
    catch(const SettingNotFoundException &nfex)
    {
      input_info_->nwrite = 10;
    }
    try
    {
      bool flag = cfg.lookup("diagnostics.write_field_timeseries");
      input_info_->write_field_timeseries = flag; // cfg.lookup("diagnostics.write_field_timeseries");
    }
    catch(const SettingNotFoundException &nfex)
    {
      input_info_->write_field_timeseries = TRUE;
    }
    try
    {
      input_info_->write_E = cfg.lookup("diagnostics.write_E");
    }
    catch(const SettingNotFoundException &nfex)
    {
      input_info_->write_E = TRUE;
    }
    try
    {
      input_info_->write_B = cfg.lookup("diagnostics.write_B");
    }
    catch(const SettingNotFoundException &nfex)
    {
      input_info_->write_B = TRUE;
    }
    try
    {
      input_info_->write_J = cfg.lookup("diagnostics.write_J");
    }
    catch(const SettingNotFoundException &nfex)
    {
      input_info_->write_J = TRUE;
    }
    try
    {
      input_info_->write_rho = cfg.lookup("diagnostics.write_rho");
    }
    catch(const SettingNotFoundException &nfex)
    {
      input_info_->write_rho = TRUE;
    }
    try
    {
      input_info_->write_all_fields = cfg.lookup("diagnostics.write_all_fields");
    }
    catch(const SettingNotFoundException &nfex)
    {
      input_info_->write_all_fields = FALSE;
    }
    

// dummy code to be replaced //////////////////////////////////
    input_info_->temp[0] = 1.5; 
    input_info_->temp[1] = 0.1; 

    sprintf(input_info_->distname,"distribution.dat");

    // MPI can only Bcast C strings
    char (*parts_bound)[NCHAR] = input_info_->parts_bound;
//    strcpy(parts_bound[0], "reflecting"); // x -> Left
//    strcpy(parts_bound[1], "reflecting"); // x -> Right
    strcpy(parts_bound[0], "periodic"); // x -> Left
    strcpy(parts_bound[1], "periodic"); // x -> Right
    strcpy(parts_bound[2], "periodic"); // y -> Left
    strcpy(parts_bound[3], "periodic"); // y -> Right
    strcpy(parts_bound[4], "periodic"); // z -> Left
    strcpy(parts_bound[5], "periodic"); // z -> Right
 
    char (*fields_bound)[NCHAR] = input_info_->fields_bound;
//    strcpy(fields_bound[0], "conducting"); // x -> Left
//    strcpy(fields_bound[1], "conducting"); // x -> Right
    strcpy(fields_bound[0], "periodic"); // x -> Left
    strcpy(fields_bound[1], "periodic"); // x -> Right
    strcpy(fields_bound[2], "periodic"); // y -> Left
    strcpy(fields_bound[3], "periodic"); // y -> Right
    strcpy(fields_bound[4], "periodic"); // z -> Left
    strcpy(fields_bound[5], "periodic"); // z -> Right
////////////////////////////////////////////////////////////////

    return 0;
}

