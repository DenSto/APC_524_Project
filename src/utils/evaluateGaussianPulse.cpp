#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "GaussianPulse.hpp"

//! compute superposition of gaussian pulses. t: current time. dt: phase shift 
double GaussianPulses(double t, double dt, Gaussian_Pulses_t *pulses){
    int nwaves = pulses->nwaves;
    double *amp = pulses->peakamps;
    double *pha = pulses->phases;
    double *omg = pulses->omegas;
    double *del = pulses->delays;
    double *dwt  = pulses->invWidths;
    
    assert(amp!=NULL);
    assert(pha!=NULL);
    assert(omg!=NULL);
    assert(del!=NULL);
    assert(dwt!=NULL);

    double val = 0.0;
    double tmp;
    for(int i=0;i<nwaves;i++){
       tmp = amp[i]*cos(omg[i]*(t+dt)+pha[i]*M_PI/180.0);
       val+= tmp*exp(-pow((t-del[i])*dwt[i],2));
    }

    return val;

};
