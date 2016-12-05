#ifndef PARTICLE_UTILS_HPP
#define PARTICLE_UTILS_HPP

#include "particle.hpp"
#include "../grid/grid.hpp"

// Function to use are a comparator in std::sort for std::vec<Particle>
// Idea: Particle is sorted by outer index first 
// (i.e. particle at [2][12][43] should be closer to beginning of array than
// particle at [24][12][43])
class Particle_Compare {
	double x0_;
	double y0_;
	double z0_;
	double idx_;
	double idy_;
	double idz_;
	public:
	Particle_Compare(Grid* grid) : 
		x0_(grid->getx0()),
		y0_(grid->gety0()),
		z0_(grid->getz0()),

		idx_(1.0/grid->getdx()),
		idy_(1.0/grid->getdy()),
		idz_(1.0/grid->getdz())
	{
	}

	bool operator()(Particle const* a, Particle const* b) const {
		int ai,bi,aj,bj,ak,bk;

		ai=(int)((a->x1 - x0_)*idx_);
		bi=(int)((b->x1 - x0_)*idx_);
		if (ai > bi) return 1;
		if (ai < bi) return 0;

		aj=(int)((a->x2 - y0_)*idy_);
		bj=(int)((b->x2 - y0_)*idy_);
		if (aj > bj) return 1;
		if (aj < bj) return 0;

		ak=(int)((a->x3 - z0_)*idz_);
		bk=(int)((b->x3 - z0_)*idz_);
		if (ak > bk) return 1;
		if (ak < bk) return 0;

		// if all indices are equal, return 0;
		return 0;
	}
};

#endif
