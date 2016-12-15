#ifndef PARTICLE_HANDLER_HPP
#define PARTICLE_HANDLER_HPP

//! Class that handles all particle-relevant operations.
/*!
        Particle handler handles all the particle operations. This includes deposition,
        boundary conditions, particle pushing, and communication between MPI nodes if
        needed
*/
#include <vector>
#include <stdio.h>
#include "../IO/IO.hpp"
#include "../Domain/domain.hpp"
#include "../grid/grid.hpp"
#include "particle_utils.hpp"
#include "../pusher/pusher.hpp"
#include "../boundaries/boundary_particles.hpp"
#include "interpolate.hpp"
#include "deposit.hpp"

class Particle_Handler {
public:
  Particle_Handler(long np); // list of np particles and their fields
  ~Particle_Handler();
  void Load(Input_Info_t info, Domain* domain, int restart); // Initialize particles
  void Push(double dt);   // Push all particles
  void Pass();            // Pass particles accross MPI boundary
  long nParticles();

  void incrementNParticles(int inc);

  void SortParticles(Particle_Compare comp); // quicksort particle list

  void setPusher(Pusher* pusher) {pusher_=pusher;};
  void clearGhosts();	// remove all ghost particles in particle list

  void InterpolateEB(Grid* grid);
  //void depositRhoJ(Grid *grids, double dt); // deposit current and charge density from particles to grid
  // depositRhoJ should not explicitly depends on dt, see deposit.cpp for comments 
  void depositRhoJ(Grid *grids); // deposit current and charge density from particles to grid

  std::vector<Particle> getParticleVector(){return parts_;}

  double maxVelocity(void); // return maximum velocity of particles 
                                  // to determine size of time steps 

  void setParticleBoundaries(BC_Particle** bc){boundaries_=bc;}
  void executeParticleBoundaryConditions();
private:
  long np_;
  BC_Particle** boundaries_; /* Particle Boundary Conditions */
  std::vector<Particle> parts_;    /* Vector of particles */
  Pusher* pusher_;
};

#endif
