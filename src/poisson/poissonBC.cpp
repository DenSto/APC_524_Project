#include <iostream>
#include <cstdlib>
#include <cmath>
#include <assert.h>
#include "../globals.hpp"
#include "poissonBC.hpp"

using namespace std;

PoissonBC::PoissonBC(int side, Input_Info_t *input_info){

   side_ = side;
   assert(abs(side)<=3 && abs(side)>=1);

   dim_ = abs(side_)-1; //x:0, y:1, z:2

   // convert side (-3,-2,-1,1,2,3) to index (4,2,0,1,3,5) 
   int index = (int)(2*abs((double)side_+0.25)-1.5);
   
   // load input info relevant to the boundary specified by side
   phiA_[0] = input_info->bound_phi[index];
   phiA_[1] = input_info->bound_Ax[index];
   phiA_[2] = input_info->bound_Ay[index];
   phiA_[3] = input_info->bound_Az[index];

   for(int i=0;i<4;i++){ifLoad_[i]=0;} // initialally no field is loaded
   
   if(debug>1) cerr << "rank=" << rank_MPI 
                    << ": side=" << side_
                    << ": finish poissonBC constructor."
                    << " Load from index " << index 
                    << " with values: "<<phiA_[0]<<", "<<phiA_[1]
                    << ", "<<phiA_[2]<<", "<<phiA_[3]<< endl;

};


PoissonBC::~PoissonBC(void){
}

void PoissonBC::applyBCs (double fieldID, double iternum, Grid *grids) {

   int sendID = (int)fieldID;
   int iter  = (int)iternum % 2; // iter==0, supply boundry condition to u1 fields
                                // iter!=0, supply boundary cndition to u2 fields

   // determine which field to apply boundary value to
   int which=-1; 
   if(sendID==grids->getFieldID("phi"))which=0;
   if(sendID==grids->getFieldID("Ax")) which=1;
   if(sendID==grids->getFieldID("Ay")) which=2;
   if(sendID==grids->getFieldID("Az")) which=3;
   if(debug>2)cerr<<"rank="<<rank_MPI<<": sendID="<<sendID<<", field="<<which
                  <<". phiID="<<grids->getFieldID("phi")
                  <<", AxID=" <<grids->getFieldID("Ax")
                  <<", AyID=" <<grids->getFieldID("Ay")
                  <<", AzID=" <<grids->getFieldID("Az")<<endl;
   assert(which>=0 && which <=3);

   // check previous fields have been loaded for at least once
   int ind=0;
   while(ind<which){
      if(debug>2)cerr<<"rank="<<rank_MPI<<": check previous Load: "
                     <<"which=" << which <<", ind="<<ind
                     <<", ifLoad="<<ifLoad_[ind]<<endl;
      assert(ifLoad_[ind]>0);
      ind+=1;
   }

   if(debug>2)cerr<<"rank="<<rank_MPI<<": check current Load: "
                  <<"ifLoad["<<which<<"]="<<ifLoad_[which]
                  <<", and iter="<<iter<<endl;

   // if fields have not been loaded, then load values
   if(ifLoad_[which]<2){
       if(debug>2)cerr<<"rank="<<rank_MPI<<": sendID="<<sendID
                      <<", load field="<<which
                      <<", for "<<ifLoad_[which]<<" times\n";

       // use ghost cell to load boundary values, likely to be inefficient
       size_ = grids->getGhostVecSize(sendID);
       ghost_ = new double[size_];
       for(int i=0;i<size_;i++){ghost_[i] = phiA_[which];}
       
       // set ghost vec
       grids->setGhostVec(side_,ghost_,sendID + 4*iter,0); // 0:replace 
    
       ifLoad_[which]+=1;       
       delete [] ghost_;
    }else{
       if(debug>2)cerr<<"rank="<<rank_MPI<<": sendID="<<sendID
                      <<", no need to load field="<<which<<" again!\n";
    }

};
