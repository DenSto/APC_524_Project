#ifndef GAUSSIAN_PULSES_HPP
#define GAUSSIAN_PULSES_HPP

/*! This parameter struct generating field boudary values.             
    The values are superpositions of Gaussian pulses of the form:   
      peakamps \times cos( omegas * t + phases)                      
      \times exp(-(t-delays)^2*invWidths^2)                         
    The coefficients peakamps, omegas, phases, delays,and invWidths 
    are double arrays of length nwaves. 
*/

typedef struct {
    int nwaves;
    double *peakamps;
    double *omegas;
    double *phases;
    double *delays;
    double *invWidths;
} Gaussian_Pulses_t;

double GaussianPulses(double t, double dt, Gaussian_Pulses_t *pulses);

#endif
