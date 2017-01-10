#include <iostream>
#include <cstdlib>
#include <cmath>
#include <assert.h>
#include "../globals.hpp"
#include "lightBC.hpp"

using namespace std;

LightBC::LightBC(int side, Input_Info_t *input_info){

   side_ = side;
   assert(abs(side)<=3 && abs(side)>=1);
   sign_=side_/abs(side_);
   dim_ = abs(side_)-1; //x:0, y:1, z:2

   // load background fields
   E0_ = input_info->E0;
   B0_ = input_info->B0;

   // load input info relevant to the boundary specified by side
   int nwaves_tot  = input_info->nwaves;
   int *inSide = input_info->inSide;
   int *inPolE = input_info->inPolE; 
   int *registry = new int[3*nwaves_tot];// tmp variable storing 3 polarizations

   for(int i=0;i<3;i++){nwaves_[i] = 0;};
   int ind,pid;
   for(int i=0;i<nwaves_tot;i++){
      if(inSide[i]==side){
        pid = inPolE[i]-1; //polarization x:0, y:1, z:2
        ind = nwaves_tot*pid+nwaves_[pid];
        registry[ind]=i;
        nwaves_[pid]+=1;
      }
   }
   if(debug) cerr << "rank=" << rank_MPI << ": side " << side << " is injected with: " << endl 
                    << "  "<< nwaves_[0] << " waves injected with x polarization; " << endl
                    << "  "<< nwaves_[1] << " waves injected with y polarization; " << endl
                    << "  "<< nwaves_[2] << " waves injected with z polarization.\n";

   // Load pulses in each of the 3 polarizations
   int nw;
   for(pid=0;pid<3;pid++){
      nw = nwaves_[pid];    
      if(nw>0){// there are waves with this polarization
         gaussian_pulses_[pid] = new Gaussian_Pulses_t;
         gaussian_pulses_[pid]->nwaves = nw; 
         
         // allocate arrays
         gaussian_pulses_[pid]->peakamps  = new double[nw]; 
         gaussian_pulses_[pid]->omegas    = new double[nw]; 
         gaussian_pulses_[pid]->phases    = new double[nw]; 
         gaussian_pulses_[pid]->delays    = new double[nw]; 
         gaussian_pulses_[pid]->invWidths = new double[nw]; 

         // load data      
         for(int i=0;i<nw;i++){
            ind = registry[pid*nwaves_tot+i];
            gaussian_pulses_[pid]->peakamps[i]  = input_info->peakamps[ind]; 
            gaussian_pulses_[pid]->omegas[i]    = input_info->omegas[ind]; 
            gaussian_pulses_[pid]->phases[i]    = input_info->phases[ind]; 
            gaussian_pulses_[pid]->delays[i]    = input_info->delays[ind]; 
            gaussian_pulses_[pid]->invWidths[i] = input_info->invWidths[ind]; 
         }
      // the pointers are not assigned in pol directions where there is no wave  
      }
   }
   if(debug>1) cerr << "rank=" << rank_MPI << ": finish loading pulses\n";

   // register transverse polarizations 
   eString_.assign(3,"Ei");
   bString_.assign(3,"Bi");
   if(debug>1){
      cerr << "rank=" << rank_MPI << ": finish reserving strings\n";
      //     << "eString size = " << eString_.size() << endl
      //     << "bString size = " << bString_.size() << endl;
   }

   switch (dim_) {
      case 0: // in x direction
         Epol_[0]=1;
         eString_[0]="Ey";
         Bpol_[0]=2;
         bString_[0]="Bz";

         Epol_[1]=2; 
         eString_[1]="Ez";
         Bpol_[1]=1;
         bString_[1]="By";

         // self
         Epol_[2]=0; 
         eString_[2]="Ex";
         Bpol_[2]=0; 
         bString_[2]="Bx";

         break;
      case 1: // in y direction
         Epol_[0]=2; 
         eString_[0]="Ez";
         Bpol_[0]=0; 
         bString_[0]="Bx";

         Epol_[1]=0; 
         eString_[1]="Ex";
         Bpol_[1]=2; 
         bString_[1]="Bz";

         // self
         Epol_[2]=1; 
         eString_[2]="Ey";
         Bpol_[2]=1; 
         bString_[2]="By";

         break;
      case 2: // in z direction 
         Epol_[0]=0; 
         eString_[0]="Ex";
         Bpol_[0]=1; 
         bString_[0]="By";

         Epol_[1]=1; 
         eString_[1]="Ey";
         Bpol_[1]=0; 
         bString_[1]="Bx";

         Epol_[2]=2; 
         eString_[2]="Ez";
         Bpol_[2]=2; 
         bString_[2]="Bz";

         break;
      default:
         cerr << "Invalid dimension in LightBC\n";
   }
   if(debug>1) cerr << "rank=" << rank_MPI << ": finish registrating transverse pulses\n";

   delete [] registry;
   if(debug>1) cerr << "rank=" << rank_MPI << ": finish lightBC constructor\n";

};


LightBC::~LightBC(void){
    int nw;
    for(int pid=0;pid<3;pid++){
       nw = nwaves_[pid];
       // delete allocated pointers
       if(nw>0){
          delete [] gaussian_pulses_[pid]->peakamps; 
          delete [] gaussian_pulses_[pid]->omegas; 
          delete [] gaussian_pulses_[pid]->phases; 
          delete [] gaussian_pulses_[pid]->delays; 
          delete [] gaussian_pulses_[pid]->invWidths; 
          delete gaussian_pulses_[pid];
       }   
    }  
}

//! Apply transverse light wave boundary condition to grid
/*!
   Uses setFieldAlongEdge method in grid to add field to grid.
*/
void LightBC::applyBCs (double t, double dt, Grid *grids) {

   // load background fields
   double E1[3],B1[3];
   for(int i=0;i<3;i++){
       E1[i] = E0_[i];
       B1[i] = B0_[i];
   }

   // add transverse wave fields
   double dx = grids->getStepSize(dim_);
   double bPhase = (dx + sign_*dt)/ 2; //inject from boundary into domain

   int bsign;
   double eFieldVal;
   double bFieldVal;
   int eid,bid;
   for(int i=0;i<3;i++){//2 transverse directions + 1 self
       eid = Epol_[i];
       bid = Bpol_[i];
       if(nwaves_[eid]>0 && eid!=dim_){// inject transverse waves on background
          eFieldVal = GaussianPulses(t,0.0,gaussian_pulses_[eid]);
          E1[eid]+=eFieldVal;

          // determine sign of B field, so that poyting flux into domain
          bsign = sign_*pow(-1,i+1);//pol is assigned in right-handed order
          bFieldVal = bsign*GaussianPulses(t,bPhase,gaussian_pulses_[eid]);//should be eid!
          B1[bid]+=bFieldVal;
       }
       grids->setFieldAlongEdge(eString_[i],dim_,side_>0, E1[eid]);
       grids->setFieldAlongEdge(bString_[i],dim_,side_>0, B1[bid]);
   }

};
