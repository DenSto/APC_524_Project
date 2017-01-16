#if USE_MPI

#include "../../globals.hpp"
#include "../../IO/input.hpp"
#include "../particles_boundary.hpp"
#include "../particle_bc_factory.hpp"
#include <vector>
#include <stdlib.h>
#include "assert.h"
#include "mpi.h"
#include "math.h"

#define DOUBLES_IN_PARTICLE 10

class BC_P_MPI : public BC_Particle {
	public:
		BC_P_MPI(Domain* domain, int dim_Index, short isLeft, std::string type);
		~BC_P_MPI();
		void computeParticleBCs(std::vector<Particle> *pl);
		int completeBC(std::vector<Particle> *pl);
	private:
		int particle_BC(Particle* p);
		double xMin_;
		double xMax_;
		int dim_index_;
		short isRight_;
		std::string type_;

		int sendRank_, recvRank_;
		double* lengthShiftRecv_,*lengthShiftSend_;
		long toSend_, toReceive_;
		std::vector<double> sendBuf_;
		std::vector<Particle> ghostBuf_;
		double* recvBuf_;
		int rBufSize_;
	
		Input_Info_t* info_; // keep info on particle types

		void packParticle(Particle* p);
		Particle unpackParticle(int offset);
};


BC_P_MPI::BC_P_MPI(Domain* domain, int dim_Index, short isRight, std::string type)
	:	dim_index_(dim_Index),
		isRight_(isRight),
		type_(type)
{
	assert(dim_index_ < 3);

	xMin_ = domain->getxyz0()[dim_index_];
	xMax_ = xMin_ + domain->getLxyz()[dim_index_];
	if(debug>1)fprintf(stderr,"rank=%d:dim=%d,isRight=%d,MPI_BC,xMin=%f,xMax=%f\n",
                                   rank_MPI,dim_index_,isRight_,xMin_,xMax_); 	

	toSend_ = 0;
	toReceive_ = 0;
	rBufSize_ = 1;
	recvBuf_ = (double *) malloc(sizeof(double)*rBufSize_*DOUBLES_IN_PARTICLE);
	assert(recvBuf_ != NULL);
	info_ = Part_BC_Factory::getInstance().getInfo();

	std::string periodic ("periodic");
	bool isPeriodic = (periodic.compare(type_) == 0);

	lengthShiftRecv_ = new double[3];
	lengthShiftSend_ = new double[3];
	for(int i = 0; i < 3; i++){
		lengthShiftRecv_[i] = 0.0;
		lengthShiftSend_[i] = 0.0;
	}

	int* nProc = domain->getnProcxyz();
	int* myLoc = domain->getmyijk();
	int* neigh = domain->getNeighbours();

	int partitionIndex = myLoc[dim_index_];
	short inMiddle = (partitionIndex != 0 && (partitionIndex != nProc[dim_index_]- 1));

        // handle opposite boundaries to avoid MPI deadlock
	if(inMiddle || isPeriodic){
		if(isRight_){
			sendRank_ = neigh[2*dim_index_ + 1]; //send to right
			recvRank_ = neigh[2*dim_index_ ];    //receive from left
		} else {
			sendRank_ = neigh[2*dim_index_];     //send to left
			recvRank_ = neigh[2*dim_index_ + 1]; //receive from right

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

	if(partitionIndex == 0 && isPeriodic){// Left most processor responsible for wrap-around in x
		const int* nxyz = domain->getnProcxyz();
		const double* L = domain->getLxyz();
		if(isRight_){
			lengthShiftRecv_[dim_index_] = L[dim_index_]*nxyz[dim_index_];
		} else {
			lengthShiftSend_[dim_index_] = L[dim_index_]*nxyz[dim_index_];
		}
	}
}


BC_P_MPI::~BC_P_MPI(){
	if(recvBuf_ != NULL)
		free(recvBuf_);
	delete[] lengthShiftSend_;
	delete[] lengthShiftRecv_;
}

int BC_P_MPI::completeBC(std::vector<Particle> *pl){
	// Send and receive particles.
	MPI_Request req;
	int err;

	// First send particle count. This way the receiving processor can allocate sufficient memory
	err = MPI_Isend(&toSend_, 1, MPI_LONG,sendRank_,11,MPI_COMM_WORLD,&req);
	if(err) fprintf(stderr, "rank=%d MPI_Isend error on toSend_ = %d\n",rank_MPI,err);

	err = MPI_Recv(&toReceive_, 1, MPI_LONG,recvRank_,11,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	if(err) fprintf(stderr, "rank=%d MPI_Recv error on toReceive_ = %d\n",rank_MPI,err);


	// If we're sending particles, go ahead and send them. Target processor should be waiting.
	if(toSend_ != 0){
		err=MPI_Wait(&req,MPI_STATUS_IGNORE);
		if(err) fprintf(stderr, "rank=%d (1) MPI_Wait error = %d\n",rank_MPI,err);

		err=MPI_Isend(&sendBuf_[0],toSend_*DOUBLES_IN_PARTICLE, MPI_DOUBLE, sendRank_,12,MPI_COMM_WORLD,&req);
		if(err) fprintf(stderr, "rank=%d MPI_Isend error on sendBuf_ = %d\n",rank_MPI,err);

		sendBuf_.clear();

		if(debug>1)fprintf(stderr,"rank=%d:dim=%d: Sending %ld particles.\n",
                                   rank_MPI,dim_index_,toSend_); 	
	}

	// If we're to receive particles, allocate sufficient memory and wait for target to send them.
	if(toReceive_ != 0){

		if(toReceive_ > rBufSize_){
			recvBuf_ = (double*) realloc(recvBuf_,sizeof(double)*DOUBLES_IN_PARTICLE*toReceive_);
			rBufSize_ = toReceive_;
			assert(recvBuf_ != NULL);
		}

		err=MPI_Recv(recvBuf_,toReceive_*DOUBLES_IN_PARTICLE, MPI_DOUBLE, recvRank_,12,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		if(err) fprintf(stderr, "rank=%d MPI_Recv error on recvBuf_ = %d\n",rank_MPI,err);

		for(int i = 0; i < toReceive_; i++){
			ghostBuf_.push_back(unpackParticle(i*DOUBLES_IN_PARTICLE));
		}
		pl->insert(pl->end(),ghostBuf_.begin(),ghostBuf_.end());
		ghostBuf_.clear();
	}

	err=MPI_Wait(&req,MPI_STATUS_IGNORE);
	if(err) fprintf(stderr, "rank=%d (2) MPI_Wait error = %d\n",rank_MPI,err);

	// Return the change in particle number. Returning function should increment accordingly.
	int ret = toReceive_ - toSend_;
	toSend_=0;
	toReceive_=0;

	return ret;
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
	sendBuf_.push_back(p->x[0] + lengthShiftSend_[0]);
	sendBuf_.push_back(p->x[1] + lengthShiftSend_[1]);
	sendBuf_.push_back(p->x[2] + lengthShiftSend_[2]);
	sendBuf_.push_back(p->v[0]);
	sendBuf_.push_back(p->v[1]);
	sendBuf_.push_back(p->v[2]);
	sendBuf_.push_back(p->gamma);
	sendBuf_.push_back((double)p->my_id);
	sendBuf_.push_back((double)p->initRank);
	sendBuf_.push_back((double)p->type);
	toSend_++;
}

Particle BC_P_MPI::unpackParticle(int offset){
	Particle p = new_particle();
	int i =offset;
	p.x[0] = recvBuf_[i++] - lengthShiftRecv_[0];
	p.x[1] = recvBuf_[i++] - lengthShiftRecv_[1];
	p.x[2] = recvBuf_[i++] - lengthShiftRecv_[2];
	p.v[0] = recvBuf_[i++];
	p.v[1] = recvBuf_[i++];
	p.v[2] = recvBuf_[i++];
	p.gamma= recvBuf_[i++];
	p.my_id = (long) recvBuf_[i++];
	p.initRank = (int) recvBuf_[i++];
	p.type = (short) recvBuf_[i++];

	p.isGhost = 0;

	assert(p.x[dim_index_] >= xMin_ && p.x[dim_index_] <= xMax_);

	p.m = info_->mass_ratio[p.type];
	p.q = info_->charge_ratio[p.type];
	p.isTestParticle = info_->isTestParticle[p.type];


	return p;
}


// Registers bounary condition into BC_Factory dictionary
static RegisterParticleBoundary instance("MPI", makeBCParticle<BC_P_MPI>);
#endif
