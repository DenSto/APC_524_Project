#include "../grid/grid.hpp"
#include "../particles/particle_handler.hpp"
#include "../globals.hpp"
#include "input.hpp"
#include <string.h>

//void writeoutput(Grid *grids, Particle_Handler *parts__fields); //MPI

void checkMPI(const char filestem[],double *buffer,int len);

class OutputBoxQuantities {
	public:
		OutputBoxQuantities(const char* fname, Grid* grid, Particle_Handler* handler, Input_Info_t* info);	
		~OutputBoxQuantities();

		void setParticleHandler(Particle_Handler* handler){pHandler_ = handler;}	
		void setGrid(Grid* grid){grid_ = grid;}	

		void output(double t, long i);
	private:
		Grid* grid_;
		char filename_[50];
		Particle_Handler* pHandler_;

		int dStep_,nextStep_;
		double dT_,nextT_;
};
