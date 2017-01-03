#ifndef PARTICLES_BOUNDARY_HPP
#define PARTICLES_BOUNDARY_HPP


#include "../grid/grid.hpp"
#include "../particles/particle.hpp"
#include <vector>
#include <string>


//! Class which defines a particle boundary condition
/*!
 * Boundary conditions have two stages.
 *
 * 1st stage: Cycling through particle list and determining which
 *       	  particles need to have boundary conditions applied,
 *			  then applies them.
 *
 * 2nd stage: Perform any more auxilliary computations, including MPI
 * 			  calls, creating new ghost particles, shuffling particles 
 * 			  ETC...
 */		
class BC_Particle {

	public:
		virtual ~BC_Particle() {};
		int computeParticleBCs(std::vector<Particle> *pl); // Cycle through particles. Returns the change in the number of particles
	private:
		virtual int particle_BC(Particle* p) = 0; // Returns 1 if particle is now ghost, 0 otherwise
		virtual int completeBC(std::vector<Particle> *pl) = 0; // Complete boundary conditions.
															  // returns the change in the number of particles
		std::string type_;
};

#endif
