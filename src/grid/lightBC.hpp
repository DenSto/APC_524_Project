#ifndef LIGHT_BC_HPP
#define LIGHT_BC_HPP

#include <vector>
#include <string>
#include "gridBC.hpp"
#include "../IO/input.hpp"
#include "../utils/GaussianPulse.hpp"


//! Class for supplying light-wave boundary conditions to field grid.
/*!
   Boundary conditions are of linear superpositions of the wave and constant \n
   The wave is of the form:
     peakamp * cos( omega * t + phase) * exp(-(t-delay)^2*invWidth^2) \n
     along plane perpendicular to side 
       side = -1: x left, side = +1: x right \n
       side = -2: y left, side = +2: y right \n
       side = -3: z left, side = +3: y right \n
*/
class LightBC : public GridBC {
public:
   LightBC(int side, Input_Info_t *input_info);
   ~LightBC(void);
   // t: current time, dt: for phase shift B field
   void applyBCs(double t, double dt, Grid *grids);

private:
   int sign_; //forward(+1) or backward(-1) wave propagation

   // background fields
   double *E0_,*B0_;

   // waves with 3 polarizations
   int nwaves_[3];
   Gaussian_Pulses_t *gaussian_pulses_[3]; 

   // 2 transverse polarizations + 1 self
   int Epol_[3],Bpol_[3];
   std::vector<std::string> eString_;
   std::vector<std::string> bString_;
   
};


#endif
