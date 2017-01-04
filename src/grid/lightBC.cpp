#include <string>
#include <math.h>
#include <cstdio>

#include "lightBC.hpp"
#include "grid.hpp"


LightBC::LightBC(int dim, bool edge, double amp, double omega, double phase)
{
   dim_ = dim;
   edge_ = edge;
   amp_ = amp;
   omega_ = omega;
   phase_ = phase;
};

//! Apply light wave boundary condition to grid
/*!
   Uses setFieldAlongEdge method in grid to add field to grid.
*/
void LightBC::applyBCs ( double t, double dt, Grid &grid) {

   double dx;
   double bPhase;
   std::string eString;
   std::string bString;
   double eFieldVal;
   double bFieldVal;

   // polarize according to dimension
   switch (dim_) {
      case 0: 
         eString = "Ey";
         bString = "Bz";
         break;
      case 1: 
         eString = "Ez";
         bString = "Bx";
         break;
      case 2: 
         eString = "Ex";
         bString = "By";
         break;
      default:
         printf("Invalid dimension in LightBC::applyBCs\n");
   }

   dx = grid.getStepSize(dim_);

   bPhase = ( dx / omega_ - dt * omega_ ) / 2;

   eFieldVal = amp_ * cos( omega_ * t + phase_ );
   grid.setFieldAlongEdge( eString, dim_, edge_, eFieldVal);

   // make into "incoming" light wave regardless of edge
   if (edge_ == 0) {
      bFieldVal = amp_ * cos( omega_ * t + phase_ + bPhase );
   } else if (edge_ == 1) {
      bFieldVal = - amp_ * cos( omega_ * t + phase_ + bPhase );
   } else {
      printf("Invalid edge in LightBC::applyBCs\n");
   }

   grid.setFieldAlongEdge( bString, dim_, edge_, bFieldVal);
};
