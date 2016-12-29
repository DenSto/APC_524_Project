#include <stdio.h>
#include <assert.h>
#include <string.h>

#if USE_MPI
#include "mpi.h"
#endif

#include "input.hpp"
#include "../globals.hpp"

Input::Input(void){
    input_info_ = new Input_Info_t; 
}

Input::~Input(void){

    assert(input_info_->mass_ratio!=NULL);
    assert(input_info_->charge_ratio!=NULL);
    assert(input_info_->dens_frac!=NULL);
    assert(input_info_->temp!=NULL);

    delete [] input_info_->mass_ratio;
    delete [] input_info_->charge_ratio;
    delete [] input_info_->dens_frac;
    delete [] input_info_->temp;

    delete input_info_;
}

void Input::mallocinfo(int nspecies){

    assert(nspecies>0);
    //fprintf(stderr,"rank=%d: mallocinfo\n",rank_MPI);

    double *mass = new double[nspecies];
    double *charge = new double[nspecies];
    double *dens= new double[nspecies];
    double *temp = new double[nspecies];

    input_info_->mass_ratio = mass;
    input_info_->charge_ratio = charge;
    input_info_->dens_frac= dens;
    input_info_->temp = temp;

    assert(input_info_->mass_ratio!=NULL);
    assert(input_info_->charge_ratio!=NULL);
    assert(input_info_->dens_frac!=NULL);
    assert(input_info_->temp!=NULL);

    //fprintf(stderr,"rank=%d: finish mallocinfo, nspecies=%d\n",rank_MPI,nspecies);
    for(int i=0;i<nspecies;i++){
      fprintf(stderr,"rank=%d,mass,charge,dens[%d]=%f,%f,%f\n",
                      rank_MPI,i,mass[i],charge[i],dens[i]);
    }
}

Input_Info_t* Input::getinfo(void){
    return input_info_;
} 

