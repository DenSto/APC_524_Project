#include <iostream>
#include <cstdlib>
#include <cmath>
#include <assert.h>
#include "../globals.hpp"
#include "poissonBC.hpp"

using namespace std;

PoissonBC::PoissonBC(int side, Input_Info_t *input_info, Grid *grids){

   assert(abs(side)<=3 && abs(side)>=1);
   dim_ = abs(side)-1; //x:0, y:1, z:2

   grids_ = grids;
   int nGhosts = grids_->getnGhosts();
   int *nxyzReal = grids_->getnxyzReal();
   if(side<0){//left side
      offset_ = nGhosts; //physical
   }else{//right side
      offset_ = nGhosts + nxyzReal[dim_];//ghost
   }

   fieldPtr_ = grids_->getFieldPtr();
   
   // convert side (-3,-2,-1,1,2,3) to index (4,2,0,1,3,5) 
   int ind = (int)(2*abs((double)side+0.25)-1.5);
   
   // load input info relevant to the boundary specified by side
   phiA_[0] = input_info->bound_phi[ind];
   phiA_[1] = input_info->bound_Ax[ind];
   phiA_[2] = input_info->bound_Ay[ind];
   phiA_[3] = input_info->bound_Az[ind];
   
   if(debug>1) cerr << "rank=" << rank_MPI 
                    << ": side=" << side
                    << ": poissonBC offset ="<<offset_
                    << ", Load from index " << ind 
                    << " with values: "<<phiA_[0]<<", "<<phiA_[1]
                    << ", "<<phiA_[2]<<", "<<phiA_[3]<< endl;

};


PoissonBC::~PoissonBC(void){
}

//! conver fieldID of phi and A to index 0,1,2,3
int PoissonBC::fieldIDToIndex(int sendID){
   int fid;
   
   if(sendID==grids_->getFieldID("phi1") || sendID==grids_->getFieldID("phi2")){
      fid = 0;
   }else if(sendID==grids_->getFieldID("Ax1") || sendID==grids_->getFieldID("Ax2")){
      fid = 1;
   }else if(sendID==grids_->getFieldID("Ay1") || sendID==grids_->getFieldID("Ay2")){
      fid = 2;
   }else if(sendID==grids_->getFieldID("Az1") || sendID==grids_->getFieldID("Az2")){
      fid = 3;
   }else{
      fid = -1;
   }
   
   return fid;
}

void PoissonBC::applyBCs (int sendID) {

   // load field
   fid_ = fieldIDToIndex(sendID);
   if(debug>1)cerr<<"rank="<<rank_MPI<<": sendID="<<sendID<<", fid="<<fid_ <<endl;
   assert(fid_>=0);

   fieldVal_ = phiA_[fid_]; 
   field_ = fieldPtr_[sendID];      
 
   // set field values 
   grids_->setFieldInPlane(dim_,offset_,field_,fieldVal_);

};
