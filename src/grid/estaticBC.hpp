#ifndef ESTATIC_BC_HPP
#define ESTATIC_BC_HPP

#include "gridBC.hpp"
#include "../IO/input.hpp"
#include "../utils/GaussianPulse.hpp"

//! Class for supplying electrostatic boundary conditions in a single field to field grid.
/*!
   Boundary conditions are of linear superpositions of the wave and constant \n
   The wave is of the form:
     peakamp * cos( omega * t + phase) * exp(-(t-delay)^2*invWidth^2) \n
     along plane perpendicular to side 
       side = -1: x left, side = +1: x right \n
       side = -2: y left, side = +2: y right \n
       side = -3: z left, side = +3: y right \n
*/
class ElectroStaticBC : public GridBC {
public:
   ElectroStaticBC(int side, Input_Info_t *input_info);
   ~ElectroStaticBC(void);
   void applyBCs(double t, Grid *grids);

private:
   int nwaves_;
   Input_Info_t *input_info_;
   Gaussian_Pulses_t *gaussian_pulses_; 
};


#endif
