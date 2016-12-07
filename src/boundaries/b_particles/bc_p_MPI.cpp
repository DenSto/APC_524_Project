#include "bc_p_MPI.hpp"

BC_P_MPI::BC_P_MPI(Particle_List* pl, double xMin, double xMax, int dim_Index) 
	:	pl_(pl),
		xMin_(xMin),
		xMax_(xMax),
		dim_Index_(dim_Index),
		{}

BC_P_MPI::~BC_P_MPI(){
}

void BC_P_MPI::completeBC(){
	//Send/Receives go here!
}

void BC_P_MPI::particle_BC(double* x, double* v, double xMin, double xMax){
	// Load particles into buffer if needed
}
