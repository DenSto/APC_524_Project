#ifndef DEPOSIT_HPP
#define DEPOSIT_HPP

#include "particle.hpp"

class Depositor{
public:
  Depositor();
  ~Depositor();
  
  void deposit_particle_J(Particle *part, double* lcell, double* cellverts, double* JObj); //deposit particle on cells' J

  void deposit_particle_Rho(Particle *part, double* lcell, double* cellverts, double* RhoObj); //deposit particle on cells' Rho
};

#endif
