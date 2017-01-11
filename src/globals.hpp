#ifndef GLOBALS_HPP
#define GLOBALS_HPP

// Units, see units.pdf for detail
// unit for charge
#define UNIT_CHARGE 0.58668774
// unit for electric field 1->299.79 KV/cm
#define UNIT_EFIELD 299.792458
// unit of thermal velocity
#define UNIT_VTH 0.002232
		
enum fieldID {E_X, E_Y,E_Z, B_X, B_Y, B_Z};

#ifdef MAIN_CPP
int rank_MPI, size_MPI;
int debug; // printf debug flag
double time_phys; // current physical time in simulation
double dt_phys; // current physical time step

#else // MAIN_CPP

extern int rank_MPI, size_MPI;
extern int debug;
extern double time_phys, dt_phys;


#endif // MAIN_CPP
#endif // GLOBALS_HPP
