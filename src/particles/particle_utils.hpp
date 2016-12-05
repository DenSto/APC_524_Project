#include "particles.hpp"
#include "../grid/grid.hpp"

// Function to use are a comparator in std::sort for std::vec<Particle>
// Idea: Particle is sorted by outer index first 
// (i.e. particle at [2][12][43] should be closer to beginning of array than
// particle at [24][12][43])
class particle_Compare {
	Grid* grid_;
	public:
	particle_Compare(Grid* grid) : grid_(grid) {}
	bool operator()(Particle const* a, Particle const* b) const {
		int ai,bi,aj,bj,ak,bk;

		ai=(int)((a->x1 - grid->x0)/grid->dx);
		bi=(int)((b->x1 - grid->x0)/grid->dx);
		if (ai > bi) return 1;
		if (ai < bi) return 0;

		aj=(int)((a->x2 - grid->y0)/grid->dy);
		bj=(int)((b->x2 - grid->y0)/grid->dy);
		if (aj > bj) return 1;
		if (aj < bj) return 0;

		ak=(int)((a->x3 - grid->z0)/grid->dz);
		bk=(int)((b->x3 - grid->z0)/grid->dz);
		if (ak > bk) return 1;
		if (ak < bk) return 0;

		// if all indices are equal, return 0;
		return 0
	}
};
