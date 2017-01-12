#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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
    printf("    The volume of each domain is %6.2e cm^3\n",domain_volume);

    // determine super particle scaling
    long np = input_info_->np;
    double n0 = np/domain_volume;
    assert(n0>0);

    double dens_phys = input_info_->dens_phys;
    if(dens_phys<0){dens_phys=-dens_phys;}
    else if(dens_phys==0){dens_phys=n0;}
    input_info_->dens_phys = dens_phys;

    double Ns = dens_phys/n0; // super ratio
    printf("    Particle density %6.2e cc is used to simulate physical density %6.2e cc\n",
                n0,dens_phys);
    printf("        The super particle scaling is %f\n",Ns);  

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
        assert(dens[i]>=0);
        cden += dens[i];
    }

    for(int i=0;i<nspec;i++){
        dens[i]/=cden;
        printf("        Species %d: density %6.2f%%, mass %9.3f, charge %6.3f\n",
                        i,100.0*dens[i],mass[i],charge[i]);
    }

    // check positivity of mass and temp
    double *temp = input_info_->temp;
    for(int i=0;i<nspec;i++){
       if(temp[i]<0)temp[i]=-temp[i];
       if(mass[i]<0)mass[i]=-mass[i];
       if(mass[i]==0){
          fprintf(stderr,"Mass cannot be zero!\n");
          err += 1;
       }
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

    /* check time resolution *****************************/
    double ttmp;
    printf("    Typical time scales in this simulation are:\n"); 

    // light transit time t=L/c, in unit of picosecond
    double time_light = Lxyz[0]/nCell[0]; 
    for(int i=1;i<3;i++){
        ttmp = Lxyz[i]/nCell[i];
        if(ttmp<time_light)time_light = ttmp;
    }
    time_light *= UNIT_TIME;
    assert(time_light>0);
    printf("        Light transit unit cell takes %6.3e ps\n",time_light); 

    // plasma frequency, in unit of 1THz=1/ps
    ttmp = 0.0;
    for(int i=0;i<nspec;i++){
        ttmp += pow(charge[i],2)*dens[i]/mass[i];
    }
    //printf("sum e^2n/m=%f,dens_phys=%f\n",ttmp,dens_phys);
    double omega_p = UNIT_FPE*sqrt(ttmp*dens_phys)*1e-9; // convert KHz to THz
    printf("        Plasma frequency is %6.3e THz\n",omega_p); 

    // maximum gyro frequency, in unit of 1THz=1/ps
    double *B0 = input_info_->B0;
    double B =0.0; // background B field
    for(int i=0;i<3;i++){B += pow(B0[i],2);}
    B = sqrt(B);
    double qm = fabs(charge[0]/mass[0]);//charge to mass ratio
    for(int i=1;i<nspec;i++){
        ttmp =fabs(charge[i]/mass[i]);
        if(ttmp>qm)qm=ttmp;
    }
    //printf("qm=%f,B=%f",qm,B);
    double omega_c = UNIT_FCE*qm*B*1e-3; // convert GHz to THz
    printf("        Maximum gyro frequency is %6.3e THz\n",omega_c); 

    // maximum boundary wave frequency, in unit of THz
    double *omegas = input_info_->omegas;
    double omega_e = 0.0;
    if(nwaves>0){
       omega_e = omegas[0];
       for(int i=1;i<nwaves;i++){
           ttmp = omegas[i];
       }
       if(ttmp>omega_e) omega_e = ttmp;
    } 
    //printf("nwaves=%d,omega=%f\n",nwaves,omega_e);
    omega_e *= UNIT_FRAD*1e-3; // convert GHz to THz 
    printf("        Maximum boundary wave frequency is %6.3e THz\n",omega_e); 


    // super particle scaling of mass, charge and temperature
    for(int i=0;i<nspec;i++){
       mass[i]  *= Ns;
       charge[i]*= Ns*UNIT_CHARGE; // also set charge to proper unit
       temp[i]  *= Ns;
    }
    printf("    Mass, charge, and temperature are rescaled for super particles.\n"); 

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
   fprintf(stderr,"rank=%d,np=%ld\n",rank,input_info->np);
   fprintf(stderr,"rank=%d,dens_phys=%f\n",rank,input_info->dens_phys);
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
      fprintf(stderr,"rank=%d,mass,charge,dens[%d]=%f,%f,%f\n",
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


