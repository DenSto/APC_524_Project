#include "IO.hpp"
#include<assert.h>

void readinput(char *fname,Input_Info_t *input_info){

    // Inserted by Yuan for testing
    // The following code is to be replaced
    input_info->nx      = 4;
    input_info->np      = 8;
    input_info->nt      = 2;
    input_info->restart = 0;
    input_info->dens    = 0.2;
    input_info->temp    = 1.5;

}

#if USE_MPI
  #include "mpi.h"

  /* MPI arguments for broadcasting input *****************/
  Input_Type::Input_Type(){
    
      count_ = 2; // two types

      lens_ = new int[count_];
      assert(lens_!=NULL);
      lens_[0]=4;
      lens_[1]=2;
      
      disps_ = new MPI_Aint[count_];
      assert(disps_!=NULL);
      disps_[0]=0;
      disps_[1]=disps_[0]+4*sizeof(int);

      types_ = new MPI_Datatype[count_];
      assert(types_!=NULL);
      types_[0]=MPI_INT;
      types_[1]=MPI_DOUBLE;
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
  
     printf("rank=%d,nx=%d\n",rank,input_info->nx);
     printf("rank=%d,np=%d\n",rank,input_info->np);
     printf("rank=%d,nt=%d\n",rank,input_info->nt);
     printf("rank=%d,restart=%d\n",rank,input_info->restart);
     printf("rank=%d,dens=%f\n",rank,input_info->dens);
     printf("rank=%d,temp=%f\n",rank,input_info->temp);
  
  }


#endif
