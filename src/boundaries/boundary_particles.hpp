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
		virtual void particle_BC(double *x, double *v, double x0, double L) = 0;
		double L_;
		double x0_;
		int dim_index_;
		short isLeft_;
		Particle_List* pl_;
}


BC_Particle::computeParticleBCs() {
	for(int i = pl_.begin(); i < pl_.end(); i++)
		particle_BC(&(pl_[i]->x[dim_index_]),&(pl_[i]->v[dim_index_]),x0_,L_);
	}
}

typedef void (*PartBndrPtr)(double *x, double *v, const double x0, const double L);


void part_periodic(double *x, double *v, const double x0, const double L);
void part_reflecting(double *x, double *v, const double x0, const double L);
void part_absorbing(double *x, double *v, const double x0, const double L);

#endif
