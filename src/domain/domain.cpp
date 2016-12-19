#include <stdio.h>
#include <assert.h>
#include <algorithm>

#include "domain.hpp"
#include "../globals.hpp"

Domain::Domain(Input_Info_t *input_info)
      : size_(size_MPI),
        rank_(rank_MPI){

       //fprintf(stderr,"rank=%d: call Domain constructor\n",rank_);
       nGhosts_ = 1;

       /* partition domain ********************************/
       nxyz_ = new int[3]; // local grid cells in each direction
       assert(nxyz_!=NULL);

       n2xyz_ = new int[3]; // local grid area in each direction
       assert(n2xyz_!=NULL);

       nProcxyz_ = new int[3];// procs in each direction
       assert(nProcxyz_!=NULL);

       //fprintf(stderr,"rank=%d:Finished allocation\n",rank_);
       int* nCell = input_info->nCell; 
       int* nProc = input_info->nProc;
       //fprintf(stderr,"rank=%d,size=%d,nproc[0]=%d,%d,%d\n",rank_,size_,nProc[0],nProc[1],nProc[2]); 
       assert(size_ == nProc[0]*nProc[1]*nProc[2]);

       // assign private variables
       n3xyz_=(nCell[0]*nCell[1]*nCell[2])/size_; // local grid volume
       //fprintf(stderr,"rank=%d,n3xyz=%d\n",rank,n3xyz_);
       assert(n3xyz_>0);

       for (int i=0;i<3;i++){
           nProcxyz_[i]=nProc[i];

           nxyz_[i]=nCell[i]/nProc[i];
           assert((nCell[i] % nProc[i])==0);

           n2xyz_[i]=n3xyz_/nxyz_[i];
           assert((n3xyz_ % nxyz_[i])==0);
       }           

       /* determine location on map ***********************/
       // my location
       myijk_ =  new int[3];
	   neighbours_ = new int[6];
       assert(myijk_!=NULL);
       RankToijk(rank_, myijk_);

       int itmp;
       // x neihbors in toroidal topology
       if(myijk_[0]>0){itmp = myijk_[0]-1;}
       else{itmp = nProcxyz_[0]-1;}
       rank_xl_ = ijkToRank(itmp,myijk_[1],myijk_[2]);
	   neighbours_[0] = rank_xl_;

       if(myijk_[0]<nProcxyz_[0]-1){itmp = myijk_[0]+1;}
       else{itmp = 0;}
       rank_xr_ = ijkToRank(itmp,myijk_[1],myijk_[2]);
	   neighbours_[1] = rank_xr_;

       // y neihbors in toroidal topology
       if(myijk_[1]>0){itmp = myijk_[1]-1;}
       else{itmp = nProcxyz_[1]-1;}
       rank_yl_ = ijkToRank(myijk_[0],itmp,myijk_[2]);
	   neighbours_[2] = rank_yl_;

       if(myijk_[1]<nProcxyz_[1]-1){itmp = myijk_[1]+1;}
       else{itmp = 0;}
       rank_yr_ = ijkToRank(myijk_[0],itmp,myijk_[2]);
	   neighbours_[3] = rank_yr_;

       // z neihbors in toroidal topology
       if(myijk_[2]>0){itmp = myijk_[2]-1;}
       else{itmp = nProcxyz_[2]-1;}
       rank_zl_ = ijkToRank(myijk_[0],myijk_[1],itmp);
	   neighbours_[4] = rank_zl_;

       if(myijk_[2]<nProcxyz_[2]-1){itmp = myijk_[2]+1;}
       else{itmp = 0;}
       rank_zr_ = ijkToRank(myijk_[0],myijk_[1],itmp);
	   neighbours_[5] = rank_zr_;


       /* determine physical domain size ******************/
       const double *Lxyz = input_info->Lxyz;
       Lxyz_ = new double[3];
       assert(Lxyz_!=NULL);

       const double *xyz0 = input_info->xyz0; 
       xyz0_ = new double[3];
       assert(xyz0_!=NULL);

       for(int i=0;i<3;i++){
           Lxyz_[i] = Lxyz[i]/nProcxyz_[i];
           xyz0_[i] = xyz0[i]+Lxyz_[i]*myijk_[i];
       }

       //fprintf(stderr,"rank=%d: end Domain constructor\n",rank_);
}

Domain::~Domain(){
    //printf("rank=%d: call Domain destructor\n",rank_);
    delete[] nxyz_;
    delete[] n2xyz_;
    delete[] xyz0_;
    delete[] Lxyz_;
    delete[] nProcxyz_;
    delete[] myijk_;
    delete[] neighbours_;
}

//! return rank for assigned i,j,k
int Domain::ijkToRank(int i, int j, int k){
    assert(i>=0);
    assert(j>=0);
    assert(k>=0);

    int myrank;
    myrank = i*nProcxyz_[1]*nProcxyz_[2];
    myrank+= j*nProcxyz_[2];
    myrank+= k;

    //fprintf(stderr,"rank=%d:i=%d,j=%d,k=%d,myrank=%d\n",rank_,i,j,k,myrank); 
    assert(myrank>=0 && myrank<=size_);   
 
    return myrank; 
}; 


//! assign value to allocated myijk[3] 
void Domain::RankToijk(int myrank, int *myijk){
    assert(myijk!=NULL);

    int nyz;
    nyz = nProcxyz_[1]*nProcxyz_[2];
    myijk[0] = myrank/nyz;

    int res;
    res = myrank % nyz;
    myijk[1] = res/nProcxyz_[2];
    myijk[2] = res % nProcxyz_[2];
   
}; 

int Domain::getnGhosts(void){
    //printf("rank=%d: call getnGhosts\n",rank_);
    return nGhosts_;
}

int* Domain::getnxyz(void){
    //printf("rank=%d: call getnxyz\n",rank_);
    return nxyz_;
}

int* Domain::getn2xyz(void){
    //printf("rank=%d: call getnxyz\n",rank_);
    return n2xyz_;
}

double* Domain::getxyz0(void){
    //printf("rank=%d: call getxyz0\n",rank_);
    return xyz0_;
}

double* Domain::getLxyz(void){
    //printf("rank=%d: call getLxyz\n",rank_);
    return Lxyz_;
}

int* Domain::getmyijk(void){
    return myijk_;
}

int* Domain::getNeighbours(void){
    return neighbours_;
}

int* Domain::getnProcxyz(void){
    return nProcxyz_;
}

int Domain::getxl(void){return rank_xl_;}    
int Domain::getxr(void){return rank_xr_;}
    
int Domain::getyl(void){return rank_yl_;}    
int Domain::getyr(void){return rank_yr_;}
    
int Domain::getzl(void){return rank_zl_;}    
int Domain::getzr(void){return rank_zr_;}
    
//! Find minimum grid size
double Domain::getmindx(void){
    double dxyz[3];
    for(int i=0;i<3;i++){
      dxyz[i]=Lxyz_[i]/nxyz_[i];
    }

    return *std::min_element(dxyz,dxyz+3);
}

void checkdomain(Domain *domain){

      int rank = rank_MPI;
 
      fprintf(stderr,"rank=%d: Check domain\n",rank);
      fprintf(stderr,"rank=%d,nGhosts=%d\n",rank,domain->getnGhosts());

      int *myijk = domain->getmyijk();
      fprintf(stderr,"rank=%d: myijk=%d,%d,%d\n",rank,myijk[0],myijk[1],myijk[2]);

      int *nxyz = domain->getnxyz();
      fprintf(stderr,"rank=%d,nxyz=%d,%d,%d\n",rank,nxyz[0],nxyz[1],nxyz[2]);
  
      int *n2xyz = domain->getn2xyz();
      fprintf(stderr,"rank=%d,n2xyz=%d,%d,%d\n",rank,n2xyz[0],n2xyz[1],n2xyz[2]);
  
      double *Lxyz = domain->getLxyz();
      fprintf(stderr,"rank=%d,Lxyz=%f,%f,%f\n",rank,Lxyz[0],Lxyz[1],Lxyz[2]);

      double *xyz0 = domain->getxyz0();
      fprintf(stderr,"rank=%d,xyz0=%f,%f,%f\n",rank,xyz0[0],xyz0[1],xyz0[2]);
  
      fprintf(stderr,"rank=%d: xl=%d,xr=%d,yl=%d,yr=%d,zl=%d,zr=%d\n",rank,
              domain->getxl(),domain->getxr(),
              domain->getyl(),domain->getyr(),
              domain->getzl(),domain->getzr());

}
    
