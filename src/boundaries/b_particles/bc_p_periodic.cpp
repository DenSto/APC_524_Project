#include "../../globals.hpp"
#include "../boundary_particles.hpp"
#include "../bc_factory.hpp"
#include "../../domain/domain.hpp"
#include <vector>


class BC_P_Periodic : public BC_Particle {
	public:
		BC_P_Periodic(Domain* domain, int dim_Index, short isLeft, std::string type);
		~BC_P_Periodic();
		void computeParticleBCs(std::vector<Particle> pl);
		int completeBC(std::vector<Particle> pl);
	private:
		int particle_BC(Particle* p);
		double xMin_;
		double xMax_;
		int dim_index_;
		short isRight_;
		std::string type_;
		std::vector<Particle> ghostBuf_;
};

BC_P_Periodic::BC_P_Periodic(Domain* domain, int dim_Index, short isLeft, std::string type) 
	:	dim_index_(dim_Index),
		isRight_((isLeft+1)%2),// factory use isLeft
		type_(type)
{
	assert(dim_index_ < 3);
	xMin_ = domain->getxyz0()[dim_index_];
	xMax_ = xMin_+domain->getLxyz()[dim_index_];
	if(debug>1)fprintf(stderr,"rank=%d:dim=%d,isRight=%d,periodic_BC,xMin=%f,xMax=%f\n",
                                   rank_MPI,dim_index_,isRight_,xMin_,xMax_); 	

}

BC_P_Periodic::~BC_P_Periodic(){
}

int BC_P_Periodic::completeBC(std::vector<Particle> pl){
	pl.insert(pl.end(),ghostBuf_.begin(),ghostBuf_.end());
	return 0;
}

int BC_P_Periodic::particle_BC(Particle* p){
// Persistent particles (don't create new ones, but move the original around)
/*
	if(p->x[dim_index_] < xMin_ && isRight_)
		p->x[dim_index_] += (xMax_-xMin_);
	if(p->x[dim_index_] > xMax_ && !isRight_)
		p->x[dim_index_] -= (xMax_-xMin_);
	return 0;
*/

// Non-persistent particles (create new particles, delete the ones in ghost cells
// at the end of the time step.
	if(p->x[dim_index_] < xMin_ && !isRight_){
		Particle newP = *p;
		newP.x[dim_index_] += (xMax_-xMin_);
		ghostBuf_.push_back(newP);
		p->isGhost = 1;
		return 1;
	}

	if(p->x[dim_index_] > xMax_ && isRight_){
		Particle newP = *p;
		newP.x[dim_index_] -= (xMax_-xMin_);
		ghostBuf_.push_back(newP);
		p->isGhost = 1;
		return 1;
	}
		
	return 0;
}

static RegisterParticleBoundary instance("periodic", makeBCParticle<BC_P_Periodic>);
