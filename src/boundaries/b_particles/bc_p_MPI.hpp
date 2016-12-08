#ifndef BC_P_MPI
#define BC_P_MPI

//! Class representing an MPI transfer condition for particles
/*!


*/


#include "../boudary_particles.hpp"
#include <vector>

class BC_P_MPI : public BC_Particle {
	public:
		BC_P_MPI(Particle_List* pl, double xMin, double xMax, int targetRank, double lengthShift);
		~BC_P_MPI();
		void computeParticleBCs();
		void completeBC();
	private:
		int particle_BC(double *x, double *v, double xMin, double xMax);
		double xMin_;
		double xMax_;
		int dim_index_;
		int targetRank_;
		double lengthShift_;
		short isLeft_;
		Particle_List* pl_;
}

#endif