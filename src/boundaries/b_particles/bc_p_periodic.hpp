#ifndef BC_P_PERIODIC
#define BC_P_PERIODIC

//! Class representing a periodic boundary condition for particles
/*!


*/


#include "../boundary_particles.hpp"

class BC_P_Periodic : public BC_Particle {
	public:
		BC_P_Periodic(Particle_List* pl, double xMin, double xMax, int dim_Index);
		~BC_P_Periodic();
		void computeParticleBCs();
		void completeBC();
	private:
		int particle_BC(double *x, double *v, double xMin, double xMax);
		double xMin_;
		double xMax_;
		int dim_index_;
		short isLeft_;
		Particle_List* pl_;
}

#endif
