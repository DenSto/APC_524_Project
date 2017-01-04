#if USE_MPI

#include "../../globals.hpp"
#include "../fields_boundary.hpp"
#include "../field_bc_factory.hpp"
#include <vector>
#include <stdlib.h>
#include "assert.h"
#include "mpi.h"

#define DOUBLES_IN_PARTICLE 10

class BC_F_MPI : public BC_Field {
	public:
		BC_F_MPI(Domain* domain, Grid *grids, int side);
		~BC_F_MPI();
		int completeBC();
//		void computeParticleBCs(std::vector<Particle> *pl);
	private:
//		int particle_BC(Particle* p);

/*		int sendRank_, recvRank_;
		double* lengthShift_;
		long toSend_, toReceive_;
		std::vector<double> sendBuf_;
		std::vector<Particle> ghostBuf_;
		double* recvBuf_;
		int rBufSize_;
		void packParticle(Particle* p);
		Particle unpackParticle(int offset);
*/
};


BC_F_MPI::BC_F_MPI(Domain* domain, Grid *grids, int side){
/*	assert(dim_index_ < 3);

	xMin_ = domain->getxyz0()[dim_index_];
	xMax_ = xMin_ + domain->getLxyz()[dim_index_];
	if(debug>1)fprintf(stderr,"rank=%d:dim=%d,isRight=%d,MPI_BC,xMin=%f,xMax=%f\n",
                                   rank_MPI,dim_index_,isRight_,xMin_,xMax_); 	

	toSend_ = 0;
	toReceive_ = 0;
	rBufSize_ = 1;
	recvBuf_ = (double *) malloc(sizeof(double)*rBufSize_*DOUBLES_IN_PARTICLE);
	assert(recvBuf_ != NULL);

	std::string periodic ("periodic");
	bool isPeriodic = (periodic.compare(type) == 0);

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
			sendRank_ = neigh[2*dim_index_];     //send to left
			recvRank_ = neigh[2*dim_index_ + 1]; //receive from right

			if(!inMiddle){// Left most processor responsible for wrap-around in x
				int* nxyz = domain->getnxyz();
				double* L = domain->getLxyz();
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
*/
}

BC_F_MPI::~BC_F_MPI(){
/*
	if(recvBuf_ != NULL)
		free(recvBuf_);
	delete[] lengthShift_;
*/
}

int BC_F_MPI::completeBC(){
/*	// Send and receive particles.
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
*/
    return 0;
}
/*
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
	sendBuf_.push_back(p->gamma);
	sendBuf_.push_back(p->q);
	sendBuf_.push_back(p->m);
	sendBuf_.push_back((double)p->my_id);
	toSend_++;
}

Particle BC_P_MPI::unpackParticle(int offset){
	Particle p;
	int i =offset;
	p.x[0] = recvBuf_[i++] - lengthShift_[0];
	p.x[1] = recvBuf_[i++] - lengthShift_[1];
	p.x[2] = recvBuf_[i++] - lengthShift_[2];
	p.v[0] = recvBuf_[i++];
	p.v[1] = recvBuf_[i++];
	p.v[2] = recvBuf_[i++];
	p.gamma= recvBuf_[i++];
	p.q    = recvBuf_[i++];
	p.m    = recvBuf_[i++];
	p.my_id = (long) recvBuf_[i++];
	p.isGhost = 0;

	return p;
}

*/
// Registers bounary condition into BC_Factory dictionary
static RegisterFieldBoundary instance("MPI", makeBCField<BC_F_MPI>);
#endif