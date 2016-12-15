#include "../boundary_particles.hpp"
#include "../bc_factory.hpp"
#include "../../domain/domain.hpp"
#include <vector>

class BC_P_MPI : public BC_Particle {
	public:
		BC_P_MPI(Domain* domain, int dim_Index, short isLeft, std::string type);
		~BC_P_MPI();
		void computeParticleBCs(std::vector<Particle> pl);
		void completeBC(std::vector<Particle> pl);
	private:
		int particle_BC(Particle* p);
		double xMin_;
		double xMax_;
		int dim_index_;
		short isLeft_;
		std::string type_;

		int targetRank_;
		double lengthShift_;
		std::vector<Particle> sendBuf_;
		std::vector<Particle> recvBuf_;
};


BC_P_MPI::BC_P_MPI(Domain* domain, int dim_Index, short isLeft, std::string type)
	:	dim_index_(dim_Index),
		isLeft_(isLeft),
		type_(type)
		{}

BC_P_MPI::~BC_P_MPI(){
}

void BC_P_MPI::completeBC(std::vector<Particle> pl){
	//Send/Receives go here!
}

int BC_P_MPI::particle_BC(Particle* p){
// Non-persistent particles (create new particles, delete the ones in ghost cells
// at the end of the time step.
/*	if(p->x[dim_index_] < xMin_ && isLeft_){
		Particle newP = *p;
		newP.x[dim_index_] += (xMax_-xMin_);
		ghostBuf_.push_back(newP);
		p->isGhost = 1;
		return 1;
	}

	if(p->x[dim_index_] > xMax_ && !isLeft_){
		Particle newP = *p;
		newP.x[dim_index_] -= (xMax_-xMin_);
		ghostBuf_.push_back(newP);
		p->isGhost = 1;
		return 1;
	}*/
}

static RegisterParticleBoundary instance("MPI", makeBCParticle<BC_P_MPI>);
