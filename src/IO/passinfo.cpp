#if USE_MPI

#include <stdio.h>
#include <assert.h>
#include <mpi.h>

#include "input.hpp"
#include "../globals.hpp"

//! Send input_info from master to other nodes 
/*! When Input_Info_t is modified, this method 
 *  needs to be modified correspondingly */
void Input::passinfo(void){

    // parameters specific to Input_Info_t
    int count = 4; // three types

    int nint = 2*NDIM+13*1+2*NWAVE+1*NSPEC; //2 of len NDIM 
                                            //+ 13 of len 1 
                                            //+ 2 of len NWAVE
                                            //+ 1 of len NSPEC
    int nlong = 2; // 2 of len 1
    int ndouble = 3*1+4*NSPEC+5*NWAVE+4*NDIM+4*2*NDIM; //3 of len 1 
                                                     //+ 4 of nspecies
                                                     //+ 5 of len NWAVE 
                                                     //+ 4 of len NDIM
                                                     //+ 4 of len 2*NDIM
    int nchar = (3+2*2*NDIM)*NCHAR; //filename 2 initial + 2 sets boundaries 

    // specify what are the MPI data types
    MPI_Datatype *types; // MPI type of elements in each block 
    types = new MPI_Datatype[count];

    types[0]=MPI_INT;
    types[1]=MPI_LONG;
    types[2]=MPI_DOUBLE;
    types[3]=MPI_CHAR;

    // specify how many variable of each type
    int *lens;  // len[i] is number of elements in each block of types[i]
    lens = new int[count];
    
    lens[0]=nint;
    lens[1]=nlong;
    lens[2]=ndouble;
    lens[3]=nchar;
 
    // specify displacements   
    MPI_Aint *disps; // disp[i] is byte displacement of each block
    disps = new MPI_Aint[count];
    
    disps[0]=0;
    disps[1]=disps[0]+nint*sizeof(int);
    disps[2]=disps[1]+nlong*sizeof(long);
    disps[3]=disps[2]+ndouble*sizeof(double);

    // commit new MPI data type
    MPI_Datatype infotype; // new type
    MPI_Type_create_struct(count,lens,disps,types,&infotype);
    MPI_Type_commit(&infotype);

    // Broadcast input_info_
    fprintf(stderr,"rank=%d:before MPI_Bcast\n",rank_MPI);
    MPI_Bcast(input_info_,1,infotype,0,MPI_COMM_WORLD);
    //MPI_Barrier(MPI_COMM_WORLD);
    fprintf(stderr,"rank=%d:after MPI_Bcast\n",rank_MPI);
    
    // Delete temporary arrays
    delete[] types;
    delete[] lens;
    delete[] disps;

}
#endif
