#ifndef INTERPOLATE_HPP
#define INTERPOLATE_HPP

#include "particle.hpp"

class Interpolator{
 public:
  Interpolator();
  ~Interpolator();
  void interpolate_fields(double* pos, double* lcell, double* cellvars, Field_part* field); //interpolate fields at particle position
};

#endif
