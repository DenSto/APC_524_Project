#include <stdio.h>
#include<assert.h>
#include "IO.hpp"

void readinput(char *fname,Input_Info_t *input_info){

    // Inserted by Yuan for testing
    // The following code is to be replaced
    input_info->nCell[0]      = 4;
    input_info->nCell[1]      = 4;
    input_info->nCell[2]      = 4;
    input_info->nt      = 2;
    input_info->nt      = 2;
    input_info->nt      = 2;
    input_info->restart = 0;

    input_info->np      = 8;

    input_info->t0      = 0.0;
    input_info->dens    = 0.2;
    input_info->temp    = 1.5;
#if USE_MPI
    input_info->nProc[0]      = 1;
    input_info->nProc[1]      = 1;
    input_info->nProc[2]      = 1;
#endif

	input_info->boundaries_particles[0] = ("periodic"); // x -> Left
	input_info->boundaries_particles[1] = ("periodic"); // x -> Right
	input_info->boundaries_particles[2] = ("periodic"); // y -> Left
	input_info->boundaries_particles[3] = ("periodic"); // y -> Right
	input_info->boundaries_particles[4] = ("periodic"); // z -> Left
	input_info->boundaries_particles[5] = ("periodic"); // z -> Right
 
    sprintf(input_info->distname,"distribution.dat");
}

void writeoutput(double t, int rank, Grid *grids, Particle_Handler *parts_fields){
    printf("rank %d: writing output files...\n",rank);
}


#if USE_MPI
  #include "mpi.h"

  /*! construct MPI arguments for broadcasting input_info. 
   *  Need to be modified if Input_Info_t is modified */
  Input_Type::Input_Type(){
    
      count_ = 4; // three types
      int nint = 3;
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
