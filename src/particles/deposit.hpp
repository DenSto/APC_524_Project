#ifndef DEPOSIT_HPP
#define DEPOSIT_HPP

#include "particle.hpp"

class Depositor{
 public:
  Depositor();
  ~Depositor();
  //void deposit_particle_RhoJ(long* cellID, double pos[][3], double* dpos, double* lcell, double cellverts[][3], double dt, double q, double RhoJObj[][12]); //deposit particle on cells' (Rho, J)
  //YShi: see comment in deposit.cpp
  void deposit_particle_RhoJ(long* cellID, Particle *part, double* lcell, double cellverts[][3], double RhoJObj[][12]); //deposit particle on cells' (Rho, J)
};

#endif
