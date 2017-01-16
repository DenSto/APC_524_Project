#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "resolve.hpp"
#include "../globals.hpp"

Resolution::Resolution(Input_Info_t *input_info){

    resolve_ = new Resolve_t;
    if(rank_MPI==0)printf("    Determining necessary temporal/spatial resolutions...\n"); 

    /* Determin time resolution *****************************/
    if(rank_MPI==0)printf("        Typical time scales in this simulation are:\n"); 
    // light transit time t=L/c, in unit of picosecond
    int *nCell   = input_info->nCell;
    double *Lxyz = input_info->Lxyz;  

    double time_light = Lxyz[0]/nCell[0]; 
    double ttmp;
    for(int i=1;i<3;i++){
        ttmp = Lxyz[i]/nCell[i];
        if(ttmp<time_light)time_light = ttmp;
    }
    time_light *= UNIT_TIME;
    assert(time_light>0);
    if(rank_MPI==0)printf("            Light transit unit cell takes %6.3e ps\n",time_light); 
    // MPI gather
    (resolve_->time).time_light = time_light;

    // Load data to determine plasma frequency scales
    double super_ratio = input_info->super_ratio;
    double dens_phys   = input_info->dens_phys;
    double n0 = dens_phys/super_ratio; // total density in simulation

    int nspec      = input_info->nspecies;
    double *mass   = input_info->mass_ratio; // mass in simulation, rescaled by super_ratio
    double *charge = input_info->charge_ratio; // charge in simulation, rescaled by super_ratio
    double *dens   = input_info->dens_frac; // normalized density fraction

    // plasma frequency, in unit of 1THz=1/ps
    ttmp = 0.0;
    for(int i=0;i<nspec;i++){
        ttmp += pow(charge[i],2)*dens[i]/mass[i];
    }
    //if(debug)printf("sum e^2n/m=%f,n0=%f,dens_phys=%f\n",ttmp,n0,dens_phys);
    double omega_p = UNIT_FPE*sqrt(ttmp*n0)*1e-9; // convert KHz to THz
    // MPI gather
    if(rank_MPI==0)printf("            Plasma frequency is %6.3e THz\n",omega_p); 

/*
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
*/

}

Resolution::~Resolution(void){
    delete resolve_;
}
