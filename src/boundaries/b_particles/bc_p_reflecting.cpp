#include "bc_p_reflecting.hpp"

BC_P_Reflecting::BC_P_Reflecting(Particle_List* pl, double xMin, double xMax, int dim_Index) 
	:	pl_(pl),
		xMin_(xMin),
		xMax_(xMax),
		dim_Index_(dim_Index),
		{}

BC_P_Reflecting::~BC_P_Reflecting(){
}

void BC_P_Reflecting::completeBC(){}

void BC_P_Reflecting::particle_BC(double* x, double* v, double xMin, double xMax){
	if(*x > xMax){
		*x = 2.0*xMax - *x;
		*v=-*v;
	} 
	if(*x < xMin){
		*x = 2.0*xMin - *x;
		*v=-*v;
	}
}
