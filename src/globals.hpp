#ifndef GLOBALS_HPP
#define GLOBALS_HPP
		
#ifdef MAIN_CPP
int rank_MPI, size_MPI;
int debug; // printf debug flag
double time_phys; // current physical time in simulation

#else // MAIN_CPP

extern int rank_MPI, size_MPI;
extern int debug;
extern double time_phys;

#endif // MAIN_CPP
#endif // GLOBALS_HPP
