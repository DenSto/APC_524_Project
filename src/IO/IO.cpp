#include <stdio.h>
#include<assert.h>
#include "IO.hpp"
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <libconfig.h++>

using namespace std;
using namespace libconfig;


void readinput(char *fname,Input_Info_t *input_info){

    Config cfg;

    // Inserted by Yuan for testing
    // The following code is to be replaced
    input_info->nCell[0]      = 8;
    input_info->nCell[1]      = 4;
    input_info->nCell[2]      = 4;

#if USE_MPI    
    input_info->nProc[0]      = 4;
    input_info->nProc[1]      = 2;
    input_info->nProc[2]      = 2;
#else
    input_info->nProc[0]      = 1;
    input_info->nProc[1]      = 1;
    input_info->nProc[2]      = 1;
#endif

    input_info->nt      = 2;
    input_info->np      = 8;

    input_info->restart = 0;

    input_info->t0      = 0.0;
    input_info->dens    = 0.2;
    input_info->temp    = 1.5;

    input_info->xyz0[0]=  0.0;
    input_info->xyz0[1]=  0.0;
    input_info->xyz0[2]=  0.0;

    input_info->Lxyz[0] = 1.5;
    input_info->Lxyz[1] = 2.0;
    input_info->Lxyz[2] = 2.5;

    sprintf(input_info->distname,"distribution.dat");

    char (*parts_bound)[32] = input_info->parts_bound;
    strcpy(parts_bound[0], "periodic"); // x -> Left
    strcpy(parts_bound[1], "periodic"); // x -> Right
    strcpy(parts_bound[2], "periodic"); // y -> Left
    strcpy(parts_bound[3], "periodic"); // y -> Right
    strcpy(parts_bound[4], "periodic"); // z -> Left
    strcpy(parts_bound[5], "periodic"); // z -> Right
 
    char (*fields_bound)[32] = input_info->fields_bound;
    strcpy(fields_bound[0], "periodic"); // x -> Left
    strcpy(fields_bound[1], "periodic"); // x -> Right
    strcpy(fields_bound[2], "periodic"); // y -> Left
    strcpy(fields_bound[3], "periodic"); // y -> Right
    strcpy(fields_bound[4], "periodic"); // z -> Left
    strcpy(fields_bound[5], "periodic"); // z -> Right
}

void writeoutput(double t, int rank, Grid *grids, Particle_Handler *parts_fields){
    //printf("rank %d: writing output files...\n",rank);
}


#if USE_MPI
  #include "mpi.h"

  /*! construct MPI arguments for broadcasting input_info. 
   *  Need to be modified if Input_Info_t is modified */
  Input_Type::Input_Type(){
    
      count_ = 4; // three types
      int nint = 2*3+2*1; //2 of len 3 + 2 of len 1
      int nlong = 1; // 1 of len 1
      int ndouble = 3*1+2*3; //3 of len 1 + 2 of len 3
      int nchar = 50+2*6*32; //char 50 + 2*6 boundaries each of 32

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
 
     fprintf(stderr,"rank=%d:checkinput\n",rank); 
     const int *nCell = input_info->nCell;
     fprintf(stderr,"rank=%d,nCell=%d,%d,%d\n",rank,nCell[0],nCell[1],nCell[2]);
 
     const int *nProc = input_info->nProc;
     fprintf(stderr,"rank=%d,nProc=%d,%d,%d\n",rank,nProc[0],nProc[1],nProc[2]);

     fprintf(stderr,"rank=%d,nt=%d\n",rank,input_info->nt);
     fprintf(stderr,"rank=%d,restart=%d\n",rank,input_info->restart);
     fprintf(stderr,"rank=%d,np=%ld\n",rank,input_info->np);

     fprintf(stderr,"rank=%d,t0=%f\n",rank,input_info->t0);
     fprintf(stderr,"rank=%d,dens=%f\n",rank,input_info->dens);
     fprintf(stderr,"rank=%d,temp=%f\n",rank,input_info->temp);

     const double *xyz0 = input_info->xyz0;
     fprintf(stderr,"rank=%d,xyz0=%f,%f,%f\n",rank,xyz0[0],xyz0[1],xyz0[2]);

     const double *Lxyz = input_info->Lxyz;
     fprintf(stderr,"rank=%d,Lxyz=%f,%f,%f\n",rank,Lxyz[0],Lxyz[1],Lxyz[2]);

     fprintf(stderr,"rank=%d,distname=%s\n",rank,input_info->distname);
     fprintf(stderr,"rank=%d,parts_bound=%s\n",rank,input_info->parts_bound[0]);
     fprintf(stderr,"rank=%d,fields_bound=%s\n",rank,input_info->fields_bound[0]);
  
  }


#endif
