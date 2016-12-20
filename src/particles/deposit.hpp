#ifndef DEPOSIT_HPP
#define DEPOSIT_HPP

#include "particle.hpp"

class Depositor{
public:
  Depositor();
  ~Depositor();
  
  void deposit_particle_J(Particle *part, double* lcell, double* cellverts, double* RhoJObj); //deposit particle on cells' J
};

#endif
