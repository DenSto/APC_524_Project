#if USE_MPI

#include "../boundary_particles.hpp"
#include "../bc_factory.hpp"
#include "../../domain/domain.hpp"
#include <vector>
#include <stdlib.h>
#include "assert.h"
#include "mpi.h"

class BC_P_MPI : public BC_Particle {
	public:
		BC_P_MPI(Domain* domain, int dim_Index, short isRight, std::string type);
		~BC_P_MPI();
		void computeParticleBCs(std::vector<Particle> pl);
		void completeBC(std::vector<Particle> pl);
	private:
		int particle_BC(Particle* p);
		double xMin_;
		double xMax_;
		int dim_index_;
		short isRight_;
		std::string type_;

		int sendRank_, recvRank_;
		double* lengthShift_;
		long toSend_, toReceive_;
		std::vector<double> sendBuf_;
		double* recvBuf_;
		void packParticle(Particle* p);
		Particle unpackParticle(int offset);
};


BC_P_MPI::BC_P_MPI(Domain* domain, int dim_Index, short isRight, std::string type)
	:	dim_index_(dim_Index),
		isRight_(isRight),
		type_(type)
		{
			assert(dim_index_ < 3);

			toSend_ = 0;
			toReceive_ = 0;

			std::string periodic ("periodic");
			bool isPeriodic = periodic.compare(type);

			lengthShift_ = new double[3];
			for(int i = 0; i < 3; i++) lengthShift_[i] = 0.0;

			int* nProc = domain->getnProcxyz();
			int* myLoc = domain->getmyijk();
			int* neigh = domain->getNeighbours();

			int partitionIndex = myLoc[dim_index_];
			short inMiddle = (partitionIndex != 0 && (partitionIndex != nProc[dim_index_]- 1));

			if(inMiddle || isPeriodic){
				if(isRight_){
					sendRank_ = neigh[2*dim_index_ + 1]; //send to right
					recvRank_ = neigh[2*dim_index_ ];    //receive from left
				} else {
					sendRank_ = neigh[2*dim_index_];     //send to right
					recvRank_ = neigh[2*dim_index_ + 1]; //receive from left

					if(!inMiddle){// Left most processor responsible for wrap-around in x
						int* nxyz = domain->getmyijk();
						int* L = domain->getmyijk();
						lengthShift_[dim_index_] = L[dim_index_]*nxyz[dim_index_];
					}
				}
			} else {
				if(isRight_){ 
					sendRank_ = neigh[2*dim_index_ + 1]; // send to right
					recvRank_ = neigh[2*dim_index_ + 1]; // receive from right
				} else {
					sendRank_ = neigh[2*dim_index_]; // send to left 
					recvRank_ = neigh[2*dim_index_]; // receive from left
				}
			}
		}


BC_P_MPI::~BC_P_MPI(){
	if(recvBuf_ != NULL)
		free(recvBuf_);
	delete[] lengthShift_;
}

void BC_P_MPI::completeBC(std::vector<Particle> pl){

	//Send/Receives go here!
	
	MPI_Isend(&toSend_, 1, MPI_LONG,sendRank_,11,MPI_COMM_WORLD,NULL);
	MPI_Recv(&toReceive_, 1, MPI_LONG,recvRank_,11,MPI_COMM_WORLD,NULL);

	if(toSend_ != 0){
		MPI_Isend(&sendBuf_[0],toSend_*MPI_P_SIZE, MPI_DOUBLE, sendRank_,12,MPI_COMM_WORLD,NULL);
	}

	if(toReceive_ != 0){
		recvBuf_ = (double*) realloc(recvBuf_,sizeof(double*)*MPI_P_SIZE*toReceive_);
		MPI_Isend(recvBuf_,toReceive_*MPI_P_SIZE, MPI_DOUBLE, recvRank_,12,MPI_COMM_WORLD,NULL);
	
		for(int i = 0; i < toReceive_; i++){
			unpackParticle(i*MPI_P_SIZE);
		}
	}

	toSend_=0;
	toReceive_=0;
}

int BC_P_MPI::particle_BC(Particle* p){
	if(p->x[dim_index_] < xMin_ && !isRight_){
		packParticle(p);
		p->isGhost = 1;
		return 1;
	}

	if(p->x[dim_index_] > xMax_ && isRight_){
		packParticle(p);
		p->isGhost = 1;
		return 1;
	}

	return 0;
}

void BC_P_MPI::packParticle(Particle* p){
	sendBuf_.push_back(p->x[0] + lengthShift_[0]);
	sendBuf_.push_back(p->x[1] + lengthShift_[1]);
	sendBuf_.push_back(p->x[2] + lengthShift_[2]);
	sendBuf_.push_back(p->v[0]);
	sendBuf_.push_back(p->v[1]);
	sendBuf_.push_back(p->v[2]);
	sendBuf_.push_back(p->q);
	sendBuf_.push_back(p->m);
	sendBuf_.push_back((double)p->my_id);
	toSend_++;
}

Particle BC_P_MPI::unpackParticle(int offset){
	Particle p;
	p.x[0] = recvBuf_[offset+0] - lengthShift_[0];
	p.x[1] = recvBuf_[offset+1] - lengthShift_[1];
	p.x[2] = recvBuf_[offset+2] - lengthShift_[2];
	p.v[0] = recvBuf_[offset+3];
	p.v[1] = recvBuf_[offset+4];
	p.v[2] = recvBuf_[offset+5];
	p.q    = recvBuf_[offset+6];
	p.m    = recvBuf_[offset+7];
	p.my_id = (long) recvBuf_[offset+8];
	p.isGhost = 0;

	return p;
}


static RegisterParticleBoundary instance("MPI", makeBCParticle<BC_P_MPI>);
#endif
