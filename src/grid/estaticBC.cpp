#include <iostream>
#include <cstdlib>
#include <cmath>
#include <assert.h>
#include "../globals.hpp"
#include "estaticBC.hpp"

using namespace std;

ElectroStaticBC::ElectroStaticBC(int side, Input_Info_t *input_info){

   side_ = side;

   input_info_ = input_info;
   assert(input_info_!=NULL);

   // load input info relevant to the boundary specified by side
   int nwaves  = input_info_->nwaves;
   int *inSide = input_info_->inSide;
   int *inPolE = input_info_->inPolE; 
   int *registry = new int[nwaves];
   int nw = 0;
   for(int i=0;i<nwaves;i++){
      if(inSide[i]==side && inPolE[i]==abs(side)){
        registry[nw]=i;
        nw+=1;
      }
   }
   nwaves_ = nw;
   if(debug>1) cerr << "rank=" << rank_MPI << ": side" << side 
                    << " has " << nwaves_ << " electrostatic waves injected\n";

   if(nwaves_>0){
      gaussian_pulses_ = new Gaussian_Pulses_t;
      gaussian_pulses_->peakamps  = new double[nwaves_]; 
      gaussian_pulses_->omegas    = new double[nwaves_]; 
      gaussian_pulses_->phases    = new double[nwaves_]; 
      gaussian_pulses_->delays    = new double[nwaves_]; 
      gaussian_pulses_->invWidths = new double[nwaves_]; 
      
      int ind;
      for(int i=0;i<nwaves_;i++){
         ind = registry[i];
         gaussian_pulses_->peakamps[i]  = input_info_->peakamps[ind]; 
         gaussian_pulses_->omegas[i]    = input_info_->omegas[ind]; 
         gaussian_pulses_->phases[i]    = input_info_->phases[ind]; 
         gaussian_pulses_->delays[i]    = input_info_->delays[ind]; 
         gaussian_pulses_->invWidths[i] = input_info_->invWidths[ind]; 
      }
   }

   delete [] registry;

};


ElectroStaticBC::~ElectroStaticBC(void){
    if(nwaves_>0){
       delete [] gaussian_pulses_->peakamps; 
       delete [] gaussian_pulses_->omegas; 
       delete [] gaussian_pulses_->phases; 
       delete [] gaussian_pulses_->delays; 
       delete [] gaussian_pulses_->invWidths; 
       
       delete gaussian_pulses_;
    }     
}

//! Inject electrostatic wave boundary condition to grid
/*!
   Uses setFieldAlongEdge method in grid to add field to grid.
*/
void ElectroStaticBC::applyBCs(double t, Grid *grids) {

   std::string str;

   // set constant background B-field
   double *B0 = input_info_->B0;
   grids->setFieldAlongEdge(str.assign("Bx"),abs(side_),side_>0, B0[0]);
   grids->setFieldAlongEdge(str.assign("By"),abs(side_),side_>0, B0[1]);
   grids->setFieldAlongEdge(str.assign("Bz"),abs(side_),side_>0, B0[2]);

   // set e field superimposed on constant background
   double fieldVal;
   double *E0 = input_info_->E0;
   double E1[3];
   for(int i=0;i<3;i++){E1[i] = E0[i];}

   if(nwaves_>0){// inject gaussian pulses on background
      fieldVal = GaussianPulses(t,gaussian_pulses_);
      E1[abs(side_)]+=fieldVal;
   }

   grids->setFieldAlongEdge(str.assign("Ex"),abs(side_),side_>0, E1[0]);
   grids->setFieldAlongEdge(str.assign("Ey"),abs(side_),side_>0, E1[1]);
   grids->setFieldAlongEdge(str.assign("Ez"),abs(side_),side_>0, E1[2]);

};
