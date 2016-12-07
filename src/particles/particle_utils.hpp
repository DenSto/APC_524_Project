#ifndef PARTICLE_UTILS_HPP
#define PARTICLE_UTILS_HPP

#include "particle.hpp"
#include "../grid/grid.hpp"

// Function to use are a comparator in std::sort for std::vec<Particle>
// Idea: Particle is sorted by outer index first 
// (i.e. particle at [2][12][43] should be closer to beginning of array than
// particle at [24][12][43])
class Particle_Compare {
	Grid* grid_;
	public:
	Particle_Compare(Grid* grid) : 
		grid_(grid)
	{
	}

	bool operator()(Particle const* a, Particle const* b) const {
		int id1 = grid_->getCellID(a->x[0],a->x[1],a->x[2]);
		int id2 = grid_->getCellID(b->x[0],b->x[1],b->x[2]);

		// if all indices are equal, return 0;
		if(id1 <= id2)
			return 0;
		else
			return 1;
	}
};

#endif
