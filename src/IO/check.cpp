#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "input.hpp"
#include "../globals.hpp"

int Input::checkinfo(void){
    return 0;
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


