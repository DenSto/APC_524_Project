#include "../grid/grid.hpp"
#include "../particles/particle_handler.hpp"
#include <string.h>

void writeoutput(double t, int rank, Grid *grids, Particle_Handler *parts__fields); //MPI

void checkMPI(const char filestem[],double *buffer,int len);
