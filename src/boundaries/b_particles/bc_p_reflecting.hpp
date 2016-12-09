#ifndef BC_P_REFLECTING
#define BC_P_REFLECTING

//! Class representing a reflecting boundary condition for particles
/*!


*/


#include "../boundary_particles.hpp"

class BC_P_Reflecting : public BC_Particle {
	public:
		BC_P_Reflecting(Particle_List* pl, double xMin, double xMax, int dim_Index);
		~BC_P_Reflecting();
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
