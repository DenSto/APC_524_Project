#ifndef PARTICLE_LIST_HPP
#define PARTICLE_LIST_HPP

#include <vector>
#include "../grid/grid.hpp"
#include "particle_utils.hpp"
#include "../pusher/pusher.hpp"

class Particle_Field_List {
    public:
        Particle_Field_List(long np); // list of np particles and their fields
        ~Particle_Field_List();
        void Load();      // Initialize particles
        void Push(double dt); // Push all particles
        void Pass();          // Pass particles accross MPI boundary
		long nParticles();

		void SortParticles(Particle_Compare comp); // quicksort particle list

		void setPusher(Pusher* pusher) {pusher_=pusher;};

                void InterpolateEB();

    private:
        long np_;
		std::vector<Particle*> parts_;    /* Vector of particles */
		Pusher* pusher_;
};

#endif
