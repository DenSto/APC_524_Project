#ifndef GLOBALS_HPP
#define GLOBALS_HPP
		
#ifdef MAIN_CPP
int rank_MPI, size_MPI;
double time_phys; // current physical time in simulation

#else // MAIN_CPP

extern int rank_MPI, size_MPI;

#endif // MAIN_CPP
#endif // GLOBALS_HPP
