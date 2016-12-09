#include <stdio.h>
#include <assert.h>
#include <algorithm>
#include "domain.hpp"

Domain::Domain(int size, int rank, Input_Info_t *input_info)
      : size_(size),
        rank_(rank){

       //printf("rank=%d: call Domain constructor\n",rank_);
       nGhosts_ = 1;

       const int* nCell = input_info->nCell; 
       nxyz_ = new int[3];
       assert(nxyz_!=NULL);
       nxyz_[0]=nCell[0];
       nxyz_[1]=nCell[1];
       nxyz_[2]=nCell[2];
#if USE_MPI
       const int* nProc = input_info->nProc; 
	   assert(size == nProc[0]*nProc[1]*nProc[2]);
       nxyz_[0]/=nProc[0];
       nxyz_[1]/=nProc[1];
       nxyz_[2]/=nProc[2];
#endif

       //const double *xyz = input_info->xyz; 
       xyz0_ = new double[3];
       assert(xyz0_!=NULL);
       xyz0_[0] = 1.0*rank;
       xyz0_[1] = 0.0;
       xyz0_[2] = 0.0;

       //const double *lxyz = input_info->lxyz;
       Lxyz_ = new double[3];
       assert(Lxyz_!=NULL);
       Lxyz_[0] = 1.0/size;
       Lxyz_[1] = 0.5*(rank+1); // for testing
       Lxyz_[2] = 0.8;

#if USE_MPI
		nProcxyz_ = new int[3];
		myLocationOnMap_ = new int[3];
        assert(nProcxyz_!=NULL);
        assert(myLocationOnMap_!=NULL);

        nProcxyz_[0] = nProc[0]; 
        nProcxyz_[1] = nProc[1]; 
        nProcxyz_[2] = nProc[2]; 

		procMap_ = (int***) malloc(nProc[0]*sizeof(int**));
		for(int i = 0; i < nProcxyz_[0]; i++){
			procMap_[i] = (int**) malloc(nProc[1]*sizeof(int*));
			for(int j = 0; j < nProcxyz_[1]; j++){
				procMap_[i][j] =(int*) malloc(nProc[2]*sizeof(int));
			}
		}

        assert(procMap_!=NULL);

		int count = 0;
		for(int i = 0; i < nProcxyz_[0]; i++){
			for(int j = 0; j < nProcxyz_[1]; j++){
				for(int k = 0; k < nProcxyz_[2]; k++){
					procMap_[i][j][k] = count;
					if(rank_ == count){
						myLocationOnMap_[0] = i;
						myLocationOnMap_[1] = j;
						myLocationOnMap_[2] = k;
					}
					count++;		
				}
			}
		}
		
#endif
  
}

Domain::~Domain(){
    //printf("rank=%d: call Domain destructor\n",rank_);
    delete[] nxyz_;
    delete[] xyz0_;
    delete[] Lxyz_;
#if USE_MPI
	delete[] nProcxyz_;
	delete[] myLocationOnMap_;
	delete[] globalnxyz_;
	delete[] globalxyz0_;
 	delete[] globalLxyz_;
	for(int i = 0; i < nProcxyz_[0]; i++){
		for(int j = 0; j < nProcxyz_[1]; j++){
			free(procMap_[i][j]);
		}
		free(procMap_[i]);
	}
	free(procMap_);
#endif
}

int Domain::getnGhosts(void){
    //printf("rank=%d: call getnGhosts\n",rank_);
    return nGhosts_;
}

int* Domain::getnxyz(void){
    //printf("rank=%d: call getnxyz\n",rank_);
    return nxyz_;
}

double* Domain::getxyz0(void){
    //printf("rank=%d: call getxyz0\n",rank_);
    return xyz0_;
}

double* Domain::getLxyz(void){
    //printf("rank=%d: call getLxyz\n",rank_);
    return Lxyz_;
}


#if USE_MPI
int* Domain::getGlobalnyxz(void){
	return globalnxyz_;
}

double* Domain::getGlobalxyz0(void){
	return globalxyz0_;
}

double* Domain::getGlobalLxyz(void){
	return globalLxyz_;
}

int* Domain::getMyLocationOnMap(){
	return myLocationOnMap_;
}

int* Domain::getnProcxyz(){
	return nProcxyz_;
}

int*** Domain::getProcMap(){
	return procMap_;
}
#endif

//! Find minimum grid size
double Domain::getmindx(void){
    double dxyz[3];
    for(int i=0;i<3;i++){
      dxyz[i]=Lxyz_[i]/nxyz_[i];
    }

    return *std::min_element(dxyz,dxyz+3);
}

void checkdomain(int rank, Domain *domain){
 
      printf("rank=%d: Check domain\n",rank);
      printf("rank=%d,nGhosts=%d\n",rank,domain->getnGhosts());

      int *nxyz = domain->getnxyz();
      printf("rank=%d,nxyz=%d,%d,%d\n",rank,nxyz[0],nxyz[1],nxyz[2]);
  
      double *xyz0 = domain->getxyz0();
      printf("rank=%d,xyz0=%f,%f,%f\n",rank,xyz0[0],xyz0[1],xyz0[2]);
  
      double *Lxyz = domain->getLxyz();
      printf("rank=%d,Lxyz=%f,%f,%f\n",rank,Lxyz[0],Lxyz[1],Lxyz[2]);
}
    
