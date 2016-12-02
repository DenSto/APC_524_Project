#ifndef IO_HPP
#define IO_HPP

typedef struct {
    int nx;
    int np;
    int nt;
    int restart;
    double density;
    double temperature;
} Input_Info_t;

// type infomation for passing Input_Info_t in MPI
#if USE_MPI
#include "mpi.h"

typedef struct {
   int count; // number of data type blocks
   int *len;  // len[i] is number of elements in each block of types[i]
   MPI_Aint *disp; // disp[i] is byte displacement of each block
   MPI_Datatype *types; // MPI type of elements in each block 
} Input_Type_t;

Input_Type_t* new_input_type();
    
#endif

void readinput(char *fname,Input_Info_t input_info);

#endif
