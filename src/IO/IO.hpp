#ifndef IO_HPP
#define IO_HPP

// NRM: why are these needed?
#include "../grid/grid.hpp"
#include "../particles/particle_handler.hpp"
#include <string.h>

//! Structure storing info in the input file
typedef struct {
    int nCell[3]; // number of cells in each direction
    int nProc[3]; // number of processors to use in each direction

    int nt; // number of time steps
    int restart; // How many previous runs?
                 // Initial run if restart = 0

    long np; // number of particles in each domain

    double t0;   // start time of simulation
    double dens; // density
    double temp; // temperature

    double xyz0[3];
    double Lxyz[3]; 

    // MPI can only send c_str!
    char distname[50]; // name of file containing distribution function 
    char parts_bound[6][32]; // particle boundary conditions for 6 sides of box
    char fields_bound[6][32];// field boundary conditions for 6 sides of the box

} Input_Info_t;

#if USE_MPI
  #include "mpi.h"

//! Type infomation for passing Input_Info_t in MPI
/*! When Input_Info_t is modified, constructor of this 
 *  class needs to be modified correspondly. */
  class Input_Type{
      public:
          Input_Type();
          ~Input_Type();

          int getcount(void);
          int *getlens(void);
          MPI_Aint *getdisps(void);
          MPI_Datatype *gettypes(void);
      
      private:
         int count_; // number of data type blocks
         int *lens_;  // len[i] is number of elements in each block of types[i]
         MPI_Aint *disps_; // disp[i] is byte displacement of each block
         MPI_Datatype *types_; // MPI type of elements in each block 

  };

#endif

int readinput(char *fname, Input_Info_t *input_info, int size);
void checkinput(int rank, Input_Info_t *input_info);

void writeoutput(double t, int rank, Grid *grids, Particle_Handler *parts__fields); //MPI

#endif
