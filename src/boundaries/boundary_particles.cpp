#include "boundary_particles.hpp"

void BC_Particle::computeParticleBCs(std::vector<Particle> pl) {
	for(std::vector<Particle>::iterator ptr = pl.begin(); ptr != pl.end(); ++ptr){
		ptr->isGhost = ptr->isGhost || 
			  particle_BC(&(*ptr));
	}
}
