#include "../../globals.hpp"
#include "../boundary_particles.hpp"
#include "../bc_factory.hpp"
#include "../../domain/domain.hpp"
#include <vector>

#define USE_GHOST


class BC_P_Periodic : public BC_Particle {
	public:
		BC_P_Periodic(Domain* domain, int dim_Index, short isRight, std::string type);
		~BC_P_Periodic();
		void computeParticleBCs(std::vector<Particle> *pl);
		int completeBC(std::vector<Particle> *pl);
	private:
		int particle_BC(Particle* p);
		double xMin_;
		double xMax_;
		int dim_index_;
		short isRight_;
		std::string type_;
#ifdef USE_GHOST
		std::vector<Particle> ghostBuf_;
#endif
};

BC_P_Periodic::BC_P_Periodic(Domain* domain, int dim_Index, short isRight, std::string type) 
	:	dim_index_(dim_Index),
		isRight_(isRight),
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

int BC_P_Periodic::completeBC(std::vector<Particle> *pl){
#ifdef USE_GHOST
	pl->insert(pl->end(),ghostBuf_.begin(),ghostBuf_.end());
	ghostBuf_.clear();
#endif
	return 0;
}

int BC_P_Periodic::particle_BC(Particle* p){
	printf("index %lf %lf %lf %d\n]",p->x[dim_index_], xMin_,xMax_,isRight_);
#ifdef USE_GHOST
// Non-persistent particles (create new particles, delete the ones in ghost cells
// at the end of the time step.
	printf("farage\n");
	if(p->x[dim_index_] < xMin_ && !isRight_){ //left boundary
		Particle newP = *p;
		newP.x[dim_index_] += (xMax_-xMin_);
		ghostBuf_.push_back(newP);
		p->isGhost = 1;
		return 1;
	}

	if(p->x[dim_index_] > xMax_ && isRight_){ // right boundary
		Particle newP = *p;
		newP.x[dim_index_] -= (xMax_-xMin_);
		ghostBuf_.push_back(newP);
		p->isGhost = 1;
		return 1;
	}
		
	return 0;
#else
// Persistent particles (don't create new ones, but move the original around)
	printf("dorf\n");
	if(p->x[dim_index_] < xMin_ && !isRight_)
		p->x[dim_index_] += (xMax_-xMin_);
	if(p->x[dim_index_] > xMax_ && isRight_)
		p->x[dim_index_] -= (xMax_-xMin_);
	return 0;
#endif
}

// Registers bounary condition into BC_Factory dictionary
static RegisterParticleBoundary instance("periodic", makeBCParticle<BC_P_Periodic>);
