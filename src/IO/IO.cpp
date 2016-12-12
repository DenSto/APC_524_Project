#include <stdio.h>
#include<assert.h>
#include "IO.hpp"
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <libconfig.h++>

using namespace std;
using namespace libconfig;


int readinput(char *fname,Input_Info_t *input_info, int size){

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
        input_info->nCell[0] = nCell[0];
        input_info->nCell[1] = nCell[1];
        input_info->nCell[2] = nCell[2];
      } else {
        input_info->nCell[0] = nCell[0];
        input_info->nCell[1] = nCell[0];
        input_info->nCell[2] = nCell[0];
        cerr << "Error: nCell is not a 3 element array in input file."
           << endl << "Assuming nCell[0]=nCell[1]=nCell[2]=" 
           << input_info->nCell[0] << "." << endl;
      }
    }
    catch(const SettingNotFoundException &nfex)
    {
      cerr << "Error: nCell not set in input file" << endl; 
      return(EXIT_FAILURE);
    }

#if USE_MPI
    try
    {
      const Setting &nProc = cfg.lookup("domain.nProc");
      if(nProc.getLength() == 3) {
        input_info->nProc[0] = nProc[0];
        input_info->nProc[1] = nProc[1];
        input_info->nProc[2] = nProc[2];
      } else {
        input_info->nProc[0] = nProc[0];
        input_info->nProc[1] = nProc[0];
        input_info->nProc[2] = nProc[0];
        cerr << "Error: nProc is not a 3 element array in input file."
           << endl << "Assuming nProc[0]=nProc[1]=nProc[2]=" 
           << input_info->nProc[0] << "." << endl;
      }
      if(input_info->nProc[0]*input_info->nProc[1]*input_info->nProc[2]!=size)
      {
        cerr << "Error: nProc layout specified in input file does not match "
             << "number of MPI processes requested." << endl;
        return(EXIT_FAILURE);
      }
    }
    catch(const SettingNotFoundException &nfex)
    {
      cerr << "Error: nCell not set in input file" << endl; 
      return(EXIT_FAILURE);
    }
#endif

    try
    {
      input_info->nt = cfg.lookup("runtime.nTimesteps");
    }
    catch(const SettingNotFoundException &nfex)
    {
      cerr << "Error: nTimesteps not set in input file" << endl; 
      return(EXIT_FAILURE);
    }
    try
    {
      input_info->t0 = cfg.lookup("runtime.startTime");
    }
    catch(const SettingNotFoundException &nfex)
    {
      cerr << "startTime not set in input file... " 
           << "Assuming startTime = 0" << endl;
      input_info->t0 = 0.;
    }
    
    try
    {
      input_info->restart = cfg.lookup("initialization.restart");
    }
    catch(const SettingNotFoundException &nfex)
    {
      cerr << "restart not set in input file... " 
           << "Assuming restart = 0" << endl;
      input_info->restart = 0.;
    }

    try
    {
      const Setting &nParticles = 
           cfg.lookup("initialization.particles.nParticles");
      if(nParticles.getType() == Setting::TypeInt) {
        int np = nParticles;
        input_info->np = (long) np;
        if(input_info->np < 0) {
          cerr << "Error: integer overflow..." 
               << "Use nParticles = #######L in input file"
               << "for long format or use scientific notation." << endl;
          return(EXIT_FAILURE); 
        }
      } else if (nParticles.getType() == Setting::TypeInt64) {
        input_info->np = nParticles;
      } else if (nParticles.getType() == Setting::TypeFloat) {
        double np = nParticles;
        input_info->np = (long) np;
      }
    }
    catch(const SettingNotFoundException &nfex)
    {
      cerr << "Error: nParticles not set in input file" << endl; 
      return(EXIT_FAILURE);
    }

    input_info->dens    = 0.2;
    input_info->temp    = 1.5;

    input_info->boundaries_particles[0] = ("periodic"); // x -> Left
    input_info->boundaries_particles[1] = ("periodic"); // x -> Right
    input_info->boundaries_particles[2] = ("periodic"); // y -> Left
    input_info->boundaries_particles[3] = ("periodic"); // y -> Right
    input_info->boundaries_particles[4] = ("periodic"); // z -> Left
    input_info->boundaries_particles[5] = ("periodic"); // z -> Right
 
    sprintf(input_info->distname,"distribution.dat");

    return 0;
}

void writeoutput(double t, int rank, Grid *grids, Particle_Handler *parts_fields){
    printf("rank %d: writing output files...\n",rank);
}


#if USE_MPI
  #include "mpi.h"

  /*! construct MPI arguments for broadcasting input_info. 
   *  Need to be modified if Input_Info_t is modified */
  Input_Type::Input_Type(){
    
      count_ = 4; // number of types
      int nint = 8;
      int nlong = 1;
      int ndouble = 3;
      int nchar = 50; //char of length 50

      // specify what are the MPI data types
      types_ = new MPI_Datatype[count_];
      assert(types_!=NULL);
      types_[0]=MPI_INT;
      types_[1]=MPI_LONG;
      types_[2]=MPI_DOUBLE;
      types_[3]=MPI_CHAR;

      // specify how many variable of each type
      lens_ = new int[count_];
      assert(lens_!=NULL);
      lens_[0]=nint;
      lens_[1]=nlong;
      lens_[2]=ndouble;
      lens_[3]=nchar;
   
      // specify displacements   
      disps_ = new MPI_Aint[count_];
      assert(disps_!=NULL);
      disps_[0]=0;
      disps_[1]=disps_[0]+nint*sizeof(int);
      disps_[2]=disps_[1]+nlong*sizeof(long);
      disps_[3]=disps_[2]+ndouble*sizeof(double);

  }

  Input_Type::~Input_Type(){
      delete[] lens_;
      delete[] disps_;
      delete[] types_;
  }

  int Input_Type::getcount(void){
      return count_;
  }

  int* Input_Type::getlens(void){
      return lens_;
  }

  MPI_Aint* Input_Type::getdisps(void){
      return disps_;
  }

  MPI_Datatype* Input_Type::gettypes(void){
      return types_;
  }

  /* Check MPI broadcast **********************************/
  void checkinput(int rank, Input_Info_t *input_info){
  
     printf("rank=%d,nx=%d\n",rank,input_info->nCell[0]);
     printf("rank=%d,np=%ld\n",rank,input_info->np);
     printf("rank=%d,nt=%d\n",rank,input_info->nt);
     printf("rank=%d,restart=%d\n",rank,input_info->restart);
     printf("rank=%d,t0=%f\n",rank,input_info->t0);
     printf("rank=%d,dens=%f\n",rank,input_info->dens);
     printf("rank=%d,temp=%f\n",rank,input_info->temp);
     printf("rank=%d,distname=%s\n",rank,input_info->distname);
  
  }


#endif
