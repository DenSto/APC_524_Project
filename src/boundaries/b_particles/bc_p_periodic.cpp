#include "bc_p_periodic.hpp"

BC_P_Periodic::BC_P_Periodic(Particle_List* pl, double xMin, double xMax, int dim_Index) 
	:	pl_(pl),
		xMin_(xMin),
		xMax_(xMax),
		dim_Index_(dim_Index),
		{}

BC_P_Periodic::~BC_P_Periodic(){
}

void BC_P_Periodic::completeBC(){}

int BC_P_Periodic::particle_BC(double* x, double* v, double xMin, double xMax){
	if(*x > xMax)
		*x-= (xMax-xMin);
	if(*x < xMin)
		*x+= (xMax-xMin);
	return 0;
}
