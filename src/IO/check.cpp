#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include "input.hpp"
#include "../globals.hpp"
#include "../domain/resolve.hpp"

//! Check input self-consistency and sufficiency, and process input information
int Input::ProcessInfo(void){
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
        assert(nCell[i]>0);
     }
     if(reset>0){
        printf("    nCell not devisible by nProc\n");
        printf("        Reset nCell to [%d, %d, %d]\n",nCell[0],nCell[1],nCell[2]);
    }

    // calculate volume of each domain
    double *Lxyz = input_info_->Lxyz;
    double total_volume = 1.0;
    for(int i=0;i<3;i++){total_volume*=Lxyz[i];}; 
    assert(total_volume>0);
    double domain_volume = total_volume/size_MPI;
    printf("    The volume of the entire simulation box is %6.2e cm^3\n",total_volume);
    printf("    The volume of each domain is %6.2e cm^3\n",domain_volume);

    /* Check run time inputs *******************************/
    // check nTimesteps
    if(input_info_->nt<0){
        input_info_->nt = 1;
        printf("    nTimesteps is negative!\n");
        printf("        Reset nTimesteps = 1\n");
    }

    /* Check particle input *******************************/
    // determine super particle scaling
    long nparticles_tot = input_info_->nparticles_tot;
    double n0 = nparticles_tot/total_volume;
    if(debug)printf("%ld %f %f\n",nparticles_tot,total_volume,n0);
    assert(n0>=0);

    double dens_phys = input_info_->dens_phys;
    if(dens_phys<0){dens_phys=-dens_phys;} // ensure positivity
    else if(dens_phys==0){dens_phys=n0;}
    input_info_->dens_phys = dens_phys;

    double super_ratio;
    if(n0>0){
        super_ratio = dens_phys/n0; 
        printf("    Particle density %6.2e cc is used to simulate physical density %6.2e cc\n",
                    n0,dens_phys);
        printf("    The super particle scaling is %f\n",super_ratio);  
    }else{
        super_ratio = -1.0;
        printf("    There is no particle in this simulation!\n");
    }
    input_info_->super_ratio = super_ratio;
 
    int nspec = input_info_->nspecies;
    printf("    There are %d species in this simulation\n",nspec);

    // check particle dens, mass
    double *dens = input_info_->dens_frac; 
    double cden = 0.0;
    for(int i=0;i<nspec;i++){
        assert(dens[i]>=0);
        cden += dens[i];
    }

    // check positivity of mass and temp
    double *mass = input_info_->mass_ratio;
    double *temp = input_info_->temp;
    for(int i=0;i<nspec;i++){
       if(temp[i]<0)temp[i]=-temp[i];
       if(mass[i]<0)mass[i]=-mass[i];
       if(mass[i]==0){
          fprintf(stderr,"Mass cannot be zero!\n");
          err += 1;
       }
    }      

    // super particle scaling of mass, charge and temperature
    if(super_ratio>0){
       double *charge = input_info_->charge_ratio;
       for(int i=0;i<nspec;i++){
          dens[i]/=cden; // normalize dens to 1 
          printf("        Species %d: normalizd density %6.2f%%. mass %9.3f, charge %6.3f\n",
                           i,100.0*dens[i],mass[i],charge[i]);
          // super particle scaling
          mass[i]  *= super_ratio;
          charge[i]*= super_ratio; 
          temp[i]  *= super_ratio;
       }
       printf("    Mass, charge, and temperature are rescaled for super particles.\n"); 
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
           fprintf(stderr,"fleft=%d,fright=%d,pleft=%d,pright=%d\n",fleft,fright,pleft,pright);
        }
    }

    // check wave injecttion
    int nwaves = input_info_->nwaves;
    assert(nwaves<NWAVE);

    // polarization specified as x:1, y:2, z:3
    int *pol = input_info_->inPolE;
    for(int i=0;i<nwaves;i++){
        assert(pol[i]>=0 && pol[i]<=3);
    }     

    return err;

}

//! Check MPI broadcast **********************************/
void checkinput(Input_Info_t *input_info){
   int rank = rank_MPI;

   /* domain **************************/
   fprintf(stderr,"rank=%d:checkinput\n",rank); 
   const int *nCell = input_info->nCell;
   fprintf(stderr,"rank=%d,nCell=%d,%d,%d\n",rank,nCell[0],nCell[1],nCell[2]);

   const int *nProc = input_info->nProc;
   fprintf(stderr,"rank=%d,nProc=%d,%d,%d\n",rank,nProc[0],nProc[1],nProc[2]);

   const double *xyz0 = input_info->xyz0;
   fprintf(stderr,"rank=%d,xyz0=%f,%f,%f\n",rank,xyz0[0],xyz0[1],xyz0[2]);

   const double *Lxyz = input_info->Lxyz;
   fprintf(stderr,"rank=%d,Lxyz=%f,%f,%f\n",rank,Lxyz[0],Lxyz[1],Lxyz[2]);

   /* runtime *************************/
   fprintf(stderr,"rank=%d,nt=%d\n",rank,input_info->nt);
   fprintf(stderr,"rank=%d,t0=%f\n",rank,input_info->t0);
   int restart = input_info->restart;
   fprintf(stderr,"rank=%d,restart=%d\n",rank,restart);

   /* particle ************************/
   fprintf(stderr,"rank=%d,nparticles_tot=%ld\n",rank,input_info->nparticles_tot);
   fprintf(stderr,"rank=%d,dens_phys=%f\n",rank,input_info->dens_phys);
   fprintf(stderr,"rank=%d,super_ratio=%f\n",rank,input_info->super_ratio);
   fprintf(stderr,"rank=%d,relativity=%d\n",rank,input_info->relativity);

   int nspecies = input_info->nspecies;
   fprintf(stderr,"rank=%d,nspecies=%d\n",rank,nspecies);

   const double *mass = input_info->mass_ratio;
   assert(mass!=NULL);
   const double *charge = input_info->charge_ratio;
   assert(charge!=NULL);
   const double *dens = input_info->dens_frac;
   assert(dens!=NULL);
   for(int i=0;i<nspecies;i++){
      fprintf(stderr,"rank=%d,rescaled mass,charge,dens[%d]=%f,%f,%f\n",
                      rank,i,mass[i],charge[i],dens[i]);
   }

   //fprintf(stderr,"rank=%d,temp=%f\n",rank,input_info->temp);

   /* fields **************************/
   fprintf(stderr,"rank=%d,electrostatic=%d\n",rank,input_info->electrostatic);
   fprintf(stderr,"rank=%d,initialization method=%s\n",rank,input_info->fields_init);
 
   const double *B0 = input_info->B0;
   fprintf(stderr,"rank=%d,B0=%f,%f,%f\n",rank,B0[0],B0[1],B0[2]);
   
   const double *E0 = input_info->E0;
   fprintf(stderr,"rank=%d,E0=%f,%f,%f\n",rank,E0[0],E0[1],E0[2]);

   int nwaves = input_info->nwaves;
   fprintf(stderr,"rank=%d,nwaves=%d\n",rank,nwaves);

   const int *inSide = input_info->inSide;
   assert(inSide!=NULL);
   const int *inPolE = input_info->inPolE;
   assert(inPolE!=NULL);

   const double *amps = input_info->peakamps;
   assert(amps!=NULL);
   const double *omegas = input_info->omegas;
   assert(omegas!=NULL);
   const double *phases = input_info->phases;
   assert(phases!=NULL);
   const double *invWidths = input_info->invWidths;
   assert(invWidths!=NULL);
   const double *delays = input_info->delays;
   assert(delays!=NULL);

   for(int i=0;i<nwaves;i++){
      fprintf(stderr,"rank=%d,side,pol[%d]=%d,%d\n",rank,i,inSide[i],inPolE[i]);
      fprintf(stderr,"rank=%d,amp,omega,phase,invWidth,delay[%d]=%f,%f,%f,%f,%f\n",
                      rank,i,amps[i],omegas[i],phases[i],invWidths[i],delays[i]);
   }
  
   // Boundary conditions for Poisson initialization
   if(restart==0 && strcmp(input_info->fields_init,"poisson")==0){
       const double *phi = input_info->bound_phi;
       assert(phi!=NULL);
       const double *Ax = input_info->bound_Ax;
       assert(phi!=NULL);
       const double *Ay = input_info->bound_Ay;
       assert(phi!=NULL);
       const double *Az = input_info->bound_Az;
       assert(phi!=NULL);
       for(int i=0;i<6;i++){
           fprintf(stderr,"rank=%d,phi,Ax,Ay,Az[%d]=%f,%f,%f,%f\n",
                           rank,i,phi[i],Ax[i],Ay[i],Az[i]);
       }
   }


   /* outputs *************************/
   // restart
   fprintf(stderr,"rank=%d,nstep_restart=%d\n",
                   rank,input_info->nstep_restart);

   // fields
   fprintf(stderr,"rank=%d,nstep_fields=%d,which_fields=%d\n",
                   rank,input_info->nstep_fields,input_info->which_fields);

   // particles
   fprintf(stderr,"rank=%d,nstep_parts=%d,output_pCount=%d\n",
                   rank,input_info->nstep_parts,input_info->output_pCount);

   /* strings *************************/
   //fprintf(stderr,"rank=%d,distname=%s\n",rank,input_info->distname);
   fprintf(stderr,"rank=%d,parts_bound=%s\n",rank,input_info->parts_bound[0]);
   fprintf(stderr,"rank=%d,fields_bound=%s\n",rank,input_info->fields_bound[0]);

}


