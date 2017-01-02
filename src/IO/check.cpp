#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "input.hpp"
#include "../globals.hpp"

//! Check input self-consistency and sufficiency
int Input::checkinfo(void){
    int err = 0;

    /* Check domain input *********************************/
    // check nCell is divisible by nProc
    int *nCell = input_info_->nCell;
    int *nProc = input_info_->nProc;
    int reset = 0;
    int resid;
    for(int i=0;i<3;i++){
        resid=nCell[i]%nProc[i];
        if(resid!=0){
           reset+=1;
           nCell[i]-=resid;
        }
     }
     if(reset>0){
        printf("    nCell not devisible by nProc\n");
        printf("        Reset nCell to [%d, %d, %d]\n",nCell[0],nCell[1],nCell[2]);
    } 
           
    /* Check run time inputs *******************************/
    // check nTimesteps
    if(input_info_->nt<0){
        input_info_->nt = 1;
        printf("    nTimesteps is negative!\n");
        printf("        Reset nTimesteps = 1\n");
    }

    /* Check particle input *******************************/
    // check particle dens, mass
    // normalize charge to proper unit 
    int nspec = input_info_->nspecies;
    printf("    There are %d species in this simulation\n",nspec);

    double *dens = input_info_->dens_frac; 
    double *mass = input_info_->mass_ratio;
    double *charge = input_info_->charge_ratio;
    double cden = 0.0;
    for(int i=0;i<nspec;i++){
        cden += dens[i];
    }

    for(int i=0;i<nspec;i++){
        dens[i]/=cden;
        printf("        Species %d: density %6.3f%%, mass %9.3f, charge %6.3f\n",
                        i,dens[i],mass[i],charge[i]);
        // normalize charge to proper unit
        charge[i]*=UNIT_CHARGE;
    }

    /* Check boundary conditions **************************/
    char (*parts_bound)[NCHAR] = input_info_->parts_bound;
    char (*fields_bound)[NCHAR] = input_info_->fields_bound;

    // check periodic boundary condition
    int fleft,fright,pleft,pright;
    for(int i=0;i<3;i++){
        fleft = strcmp(fields_bound[2*i],"periodic");            
        fright= strcmp(fields_bound[2*i+1],"periodic");            
        pleft = strcmp(parts_bound[2*i],"periodic");            
        pright= strcmp(parts_bound[2*i+1],"periodic");            
        if(fleft==0 && fright ==0 && pleft ==0 && pright ==0){
           printf("    Boundary conditions are periodic in %d-direction\n",i); 
        } else if(fleft!=0 && fright!=0 && pleft!=0 && pright!=0){
           printf("    No boundary condition is periodic in %d-direction\n",i);
        } else {
           err += 1;
           fprintf(stderr,"Inconsistent boundary conditions in %d-direction!\n",i);
        }
    }

    return err;
}

/* Check MPI broadcast **********************************/
void checkinput(Input_Info_t *input_info){
   int rank = rank_MPI;

   fprintf(stderr,"rank=%d:checkinput\n",rank); 
   const int *nCell = input_info->nCell;
   fprintf(stderr,"rank=%d,nCell=%d,%d,%d\n",rank,nCell[0],nCell[1],nCell[2]);

   const int *nProc = input_info->nProc;
   fprintf(stderr,"rank=%d,nProc=%d,%d,%d\n",rank,nProc[0],nProc[1],nProc[2]);

   fprintf(stderr,"rank=%d,nt=%d\n",rank,input_info->nt);
   fprintf(stderr,"rank=%d,restart=%d\n",rank,input_info->restart);
   fprintf(stderr,"rank=%d,np=%ld\n",rank,input_info->np);

   fprintf(stderr,"rank=%d,t0=%f\n",rank,input_info->t0);

   int nspecies = input_info->nspecies;
   fprintf(stderr,"rank=%d,nspecies=%d\n",rank,nspecies);

   const double *mass = input_info->mass_ratio;
   assert(mass!=NULL);
   const double *charge = input_info->charge_ratio;
   assert(charge!=NULL);
   const double *dens = input_info->dens_frac;
   assert(dens!=NULL);
   for(int i=0;i<nspecies;i++){
      fprintf(stderr,"rank=%d,mass,charge,dens[%d]=%f,%f,%f\n",
                      rank,i,mass[i],charge[i],dens[i]);
   }

   //fprintf(stderr,"rank=%d,temp=%f\n",rank,input_info->temp);

   const double *xyz0 = input_info->xyz0;
   fprintf(stderr,"rank=%d,xyz0=%f,%f,%f\n",rank,xyz0[0],xyz0[1],xyz0[2]);

   const double *Lxyz = input_info->Lxyz;
   fprintf(stderr,"rank=%d,Lxyz=%f,%f,%f\n",rank,Lxyz[0],Lxyz[1],Lxyz[2]);

   fprintf(stderr,"rank=%d,distname=%s\n",rank,input_info->distname);
   fprintf(stderr,"rank=%d,parts_bound=%s\n",rank,input_info->parts_bound[0]);
   fprintf(stderr,"rank=%d,fields_bound=%s\n",rank,input_info->fields_bound[0]);

}


