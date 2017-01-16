#ifndef PARTICLE_UTILS_HPP
#define PARTICLE_UTILS_HPP

#include "particle.hpp"
#include "../grid/grid.hpp"

/*! Function to use are a comparator in std::sort for std::vec<Particle>
* Idea: Particle is sorted by outer index first 
* (i.e. particle at [2][12][43] should be closer to beginning of array than
* particle at [24][12][43])
*
* At the moment, implemented very slowly!!! Should be modified two ways:
*
* 	1) Instead of comparing cell ID, compare i,j,k indice locations individually
* 	   to save time
* 	2) Bring the code to calculating i,j,k into the comparison routine
*/
class Particle_Compare {
	public:
	Particle_Compare(Grid* grid)  
	{
		idx_ = grid->getidx();
		idy_ = grid->getidy();
		idz_ = grid->getidz();

		x0_ = grid->getx0();
		y0_ = grid->gety0();
		z0_ = grid->getz0();
	}

	bool operator()(Particle const a, Particle const b) const {
// Sort by individual i,j,k
		int na,nb;	

		na = (int)((a.x[0] - x0_)*idx_);
		nb = (int)((b.x[0] - x0_)*idx_);
		if(na < nb) return 0;
		if(na > nb) return 1;

		na = (int)((a.x[1] - y0_)*idy_);
		nb = (int)((b.x[1] - y0_)*idy_);
		if(na < nb) return 0;
		if(na > nb) return 1;

		na = (int)((a.x[2] - z0_)*idz_);
		nb = (int)((b.x[2] - z0_)*idz_);
		if(na <= nb) return 0;
		if(na > nb)  return 1;

		return 0;

// Sort by cell ID (probably overkill)
/*
		int id1 = grid_->getCellID(a.x[0],a.x[1],a.x[2]);
		int id2 = grid_->getCellID(b.x[0],b.x[1],b.x[2]);

		// if all indices are equal, return 0;
		if(id1 <= id2)
			return 0;
		else
			return 1;
*/
	}

	private:
	
		double idx_,idy_,idz_;
		double x0_,y0_,z0_;
};

#endif
