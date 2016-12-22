#ifndef POISSON_HPP
#define POISSON_HPP

#include "../grid/grid.hpp"
//#include "../domain/domain.hpp"

class Poisson_Solver : protected Grid{
public:
  Poisson_Solver(int *nxyz, int nGhosts, double *xyz0, double *Lxyz);

  void initialize_poisson_fields();

protected:
  void run_poisson_solver_(double*** u0, double*** u1,double*** R,double convergenceTol,double sourceMult);

  double ***phi1_;
  double ***phi2_;

  double ***Ax1_;
  double ***Ay1_;
  double ***Az1_;
  double ***Ax2_;
  double ***Ay2_;
  double ***Az2_;
};

#endif
