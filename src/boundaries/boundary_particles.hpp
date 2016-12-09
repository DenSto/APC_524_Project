#ifndef BOUNDARY_PARTICLES_HPP
#define BOUNDARY_PARTICLES_HPP


#include "../grid/grid.hpp"
#include "../domain/domain.hpp"
#include "../particles/particle.hpp"
#include "../particles/particle_handler.hpp"
#include <vector>
#include <string>

class BC_Particle {

	public:
		virtual ~BC_Particle() {};
		void computeParticleBCs(std::vector<Particle> pl);
		virtual void completeBC() = 0;
	private:
		virtual int particle_BC(double *x, double *v, double xMin, double xMax) = 0;
		double xMin_;
		double xMax_;
		int dim_index_;
		short isLeft_;
		std::string type_;
};

class BC_P_Collection {
	public:
		BC_P_Collection();
		~BC_P_Collection();
		void executeParticleBoundaries();
	private:
		BC_Particle* boundaries_[6];
};

#endif
