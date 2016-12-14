#ifndef BOUNDARY_PARTICLES_HPP
#define BOUNDARY_PARTICLES_HPP


#include "../grid/grid.hpp"
#include "../particles/particle.hpp"
#include <vector>
#include <string>

class BC_Particle {

	public:
		virtual ~BC_Particle() {};
		void computeParticleBCs(std::vector<Particle> pl);
		virtual void completeBC(std::vector<Particle> pl) = 0;
	private:
		virtual int particle_BC(Particle* p) = 0;
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
