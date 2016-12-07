#include <stdio.h>
#include <assert.h>
#include <algorithm>
#include "domain.hpp"

Domain::Domain(int size, int rank, Input_Info_t *input_info)
      : size_(size),
        rank_(rank){

       //printf("rank=%d: call Domain constructor\n",rank_);
       nGhosts_ = 1;

       const int nx = input_info->nx; 
       nxyz_ = new int[3];
       assert(nxyz_!=NULL);
       nxyz_[0]=nx/size;
       nxyz_[1]=nx;
       nxyz_[2]=nx;

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
  
}

Domain::~Domain(){
    //printf("rank=%d: call Domain destructor\n",rank_);
    delete[] nxyz_;
    delete[] xyz0_;
    delete[] Lxyz_;
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
    
