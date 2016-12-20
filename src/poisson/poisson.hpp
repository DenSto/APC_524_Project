#ifndef POISSON_HPP
#define POISSON_HPP

#include "../domain/domain.hpp"

class Poisson_Solver{
public:
  Poisson_Solver();
  ~Poisson_Solver();

  void run_poisson_solver(Grid* grid, Domain* domain);
};

#endif
