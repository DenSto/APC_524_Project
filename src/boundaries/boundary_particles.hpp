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
		virtual int completeBC(std::vector<Particle> pl) = 0;
	private:
		virtual int particle_BC(Particle* p) = 0;
		std::string type_;
};

#endif
