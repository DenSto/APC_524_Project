#ifndef BOUNDARY_PARTICLES_HPP
#define BOUNDARY_PARTICLES_HPP


#include "../grid/grid.hpp"
#include "../domain/domain.hpp"
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
	for(auto ptr = pl_.begin(); ptr != pl_.end(); ++ptr)
		ptr->isGhost = ptr->isGhost || 
			  particle_BC(&(ptr->x[dim_index_]),&(ptr->v[dim_index_]),x0_,L_);
	}
}



class BC_P_Collection {
	public:
		BC_P_Collection();
		~BC_P_Collection();
		executeParticleBoundaries();
	private:
		BC_Particle* boundaries_[6];
}
#endif
