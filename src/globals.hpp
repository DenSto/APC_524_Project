#ifndef GLOBALS_HPP
#define GLOBALS_HPP

// unit for charge, see units.pdf for details
#define UNIT_CHARGE 0.58668774

// unit for electric field 1->299.79 KV/cm
#define UNIT_EFIELD 299.792458
		
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
