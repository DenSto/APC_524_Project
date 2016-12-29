#if USE_MPI

#include <stdio.h>
#include <assert.h>
#include <mpi.h>

#include "input.hpp"
#include "../globals.hpp"

//! Send input_info from master to other nodes 
/*! When Input_Info_t is modified, this method 
 *  needs to be modified correspondingly */
void Input::passinfo(int nspecies){

    // parameters specific to Input_Info_t
    int count = 4; // three types

    int nint = 2*3+5*1; //2 of len 3 + 5 of len 1
    int nlong = 1; // 1 of len 1
    int ndouble = 1*1+4*nspecies+2*3; //1 of len 1 + 4 of nspecies+ 2 of len 3
    int nchar = 50+2*6*32; //char 50 + 2*6 boundaries each of 32

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
    assert(nspecies==input_info_->nspecies);
    
    // Delete temporary arrays
    delete[] types;
    delete[] lens;
    delete[] disps;

}
#endif
