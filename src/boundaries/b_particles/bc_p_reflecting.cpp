#include "../../globals.hpp"
#include "../boundary_particles.hpp"
#include "../bc_factory.hpp"
#include "../../domain/domain.hpp"
#include <vector>

class BC_P_Reflecting : public BC_Particle {
	public:
		BC_P_Reflecting(Domain* domain, int dim_Index, short isLeft, std::string type);
		~BC_P_Reflecting();
		void computeParticleBCs(std::vector<Particle> pl);
		int completeBC(std::vector<Particle> pl);
	private:
		int particle_BC(Particle* p);
		double xMin_;
		double xMax_;
		int dim_index_;
		short isRight_;
		std::string type_;
};

BC_P_Reflecting::BC_P_Reflecting(Domain* domain, int dim_Index, short isLeft, std::string type) 
	:	dim_index_(dim_Index),
		isRight_((isLeft+1)%2),// factory use isLeft
		type_(type)
		{
		fprintf(stderr,"rank=%d:dim=%d,isRight=%d,reflect BC\n",rank_MPI,dim_Index,isRight_); 	

}

BC_P_Reflecting::~BC_P_Reflecting(){
}

int BC_P_Reflecting::completeBC(std::vector<Particle> pl){
	return 0;
}


int BC_P_Reflecting::particle_BC(Particle* p){
	if(p->x[dim_index_] > xMax_ && isRight_){
		p->x[dim_index_] = 2.0*xMax_ - p->x[dim_index_];
		p->v[dim_index_]=-p->v[dim_index_];
	} 
	if(p->x[dim_index_] < xMin_ && !isRight_){
		p->x[dim_index_] = 2.0*xMin_ - p->x[dim_index_];
		p->v[dim_index_]=-p->v[dim_index_];
	}
	return 0;
}

static RegisterParticleBoundary instance("reflecting", makeBCParticle<BC_P_Reflecting>);
