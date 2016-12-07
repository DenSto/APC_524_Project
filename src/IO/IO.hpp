#ifndef IO_HPP
#define IO_HPP

//! Structure storing info in the input file
typedef struct {
    int nx; // number of grids
    int nt; // number of time steps
    int restart; // How many previous runs?
                 // Initial run if restart = 0
    long np; // number of particles in each domain

    double t0;   // start time of simulation
    double dens; // density
    double temp; // temperature

    char distname[50]; // name of file containing distribution function 

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

void readinput(char *fname, Input_Info_t *input_info);

void checkinput(int rank, Input_Info_t *input_info);

#endif
