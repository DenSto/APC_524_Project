#ifndef GAUSSIAN_PULSES_HPP
#define GAUSSIAN_PULSES_HPP

/*! This is a function generating field boudary values.             
    The values are superpositions of Gaussian pulses of the form:   
      peakamps \times cos( omegas * t + phases)                      
      \times exp(-(t-delays)^2*invWidths^2)                         
    The coefficients peakamps, omegas, phases, delays,and invWidths 
    are double arrays of length nwaves. 
*/
#include <assert.h>
#include <math.h>

typedef struct {
    int nwaves;
    double *peakamps;
    double *omegas;
    double *phases;
    double *delays;
    double *invWidths;
} Gaussian_Pulses_t;

double GaussianPulses(double t, Gaussian_Pulses_t *pulses){
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
       tmp = amp[i]*cos(omg[i]*t+pha[i]*M_PI/180.0);
       val+= tmp*exp(-pow((t-del[i])*dwt[i],2));
    }

    return val;

};

#endif
