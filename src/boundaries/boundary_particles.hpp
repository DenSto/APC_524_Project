#ifndef BOUNDARY_PARTICLES_HPP
#define BOUNDARY_PARTICLES_HPP


#include "../grid/grid.hpp"
#include "../particles/particle_list.hpp"

class BC_Particle {

	public:
		virtual ~BC_Particle();
		virtual	void computeParticleBCs() = 0;
		virtual void completeBC() = 0;
	private:
		virtual int particle_BC(double *x, double *v, double xMin, double xMax) = 0;
		double xMin_;
		double xMax_;
		int dim_index_;
		short isLeft_;
		Particle_List* pl_;
}


BC_Particle::computeParticleBCs() {
	for(int i = pl_.begin(); i < pl_.end(); i++)
		pl_[i]->isGhost = pl_[i]->isGhost || 
			  particle_BC(&(pl_[i]->x[dim_index_]),&(pl_[i]->v[dim_index_]),x0_,L_);
	}
}

#endif
