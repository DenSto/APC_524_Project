#include "../boundary_particles.hpp"
#include "../bc_factory.hpp"
#include "../../domain/domain.hpp"
#include <vector>
#include "assert.h"

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
		int toSend_, toReceive_;
//		std::vector<Particle> sendBuf_;
//		std::vector<Particle> recvBuf_;
		std::vector<double> sendBuf_;
		double* recvBuf_;
		void packParticle(Particle* p);
		Particle unpackParticle(int offset);
};


BC_P_MPI::BC_P_MPI(Domain* domain, int dim_Index, short isLeft, std::string type)
	:	dim_index_(dim_Index),
		isLeft_(isLeft),
		type_(type)
		{
			assert(dim_index_ < 3);

			int* nProc = domain->getnProcxyz();
			int* myLoc = domain->getmyijk();
			int* neigh = domain->getNeighbours();
			int partitionIndex = myLoc[dim_index_];
			short inMiddle = (partitionIndex != 0 && (partitionIndex != nProc[dim_index_]- 1));
		}




BC_P_MPI::~BC_P_MPI(){
	if(recvBuf_ != NULL)
		free(recvBuf_);
}

void BC_P_MPI::completeBC(std::vector<Particle> pl){
	//Send/Receives go here!
}

int BC_P_MPI::particle_BC(Particle* p){
	if(p->x[dim_index_] < xMin_ && isLeft_){
		p->isGhost = 1;
		return 1;
	}

	if(p->x[dim_index_] > xMax_ && !isLeft_){
		p->isGhost = 1;
		return 1;
	}
}

void BC_P_MPI::packParticle(Particle* p){
	sendBuf_.push_back(p->x[0]);
	sendBuf_.push_back(p->x[1]);
	sendBuf_.push_back(p->x[2]);
	sendBuf_.push_back(p->v[0]);
	sendBuf_.push_back(p->v[1]);
	sendBuf_.push_back(p->v[2]);
	sendBuf_.push_back(p->q * p->m);
	sendBuf_.push_back((double)p->my_id);
}

Particle BC_P_MPI::unpackParticle(int offset){
	Particle p;
	p.x[0] = recvBuf_[offset+0];
	p.x[1] = recvBuf_[offset+1];
	p.x[2] = recvBuf_[offset+2];
	p.v[0] = recvBuf_[offset+3];
	p.v[1] = recvBuf_[offset+4];
	p.v[2] = recvBuf_[offset+5];
	double q,m, qm = recvBuf_[offset+6];
	if(qm < 0){
		q = -1;
		m = q*qm;
	} else {
		q = 1;
		m = qm;
	}
	p.q = q;
	p.m = m;
	p.my_id = (int) recvBuf_[offset+7];
	p.isGhost = 0;

	return p;
}

static RegisterParticleBoundary instance("MPI", makeBCParticle<BC_P_MPI>);
