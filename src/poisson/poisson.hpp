#ifndef POISSON_HPP
#define POISSON_HPP

#include "../grid/grid.hpp"
#include "../domain/domain.hpp"

//class Poisson_Solver : protected Grid{
class Poisson_Solver : public Grid{
public:
  Poisson_Solver(Domain *domain, Input_Info_t *input_info);
  ~Poisson_Solver();

  //void initialize_poisson_fields();
  void InitializeFields(Input_Info_t *input_info);

  int getGhostVecSize(const int sendID); 
  void getGhostVec(const int side, double* ghostVec, int sendID); 
  void setGhostVec(const int side, double* ghostVec, int sendID, int op); 
  void phiToE();
  void AToB();
  void constA(const double vx, const double vy, const double vz);
  void constPhi(const double v);

protected:
  void run_poisson_solver_(const int fieldID, double*** u0, double*** u1,double*** R,double convergenceTol,double sourceMult); 
  void setPoissonFieldType_();
  void setPoissonFieldPtr_();

  void phiToESingleComp_(const int fieldID, const int dir); 
  void AToBSingleComp_(const int fieldID, const int dir); 

  Domain *domain_;
  Input_Info_t *input_info_;

  double ***phi1_;
  double ***phi2_;

  double ***Ax1_;
  double ***Ay1_;
  double ***Az1_;
  double ***Ax2_;
  double ***Ay2_;
  double ***Az2_;

  const int nFieldsPoisson_; 
  
  const int phi1ID_;
  const int phi2ID_;

  const int Ax1ID_;
  const int Ay1ID_;
  const int Az1ID_;
  const int Ax2ID_;
  const int Ay2ID_;
  const int Az2ID_;

  const int xdir_; 
  const int ydir_; 
  const int zdir_; 

  // unit testing in convertFields_unittests.cc 
  friend class ConvertPrivateTest;
  FRIEND_TEST(ConvertPrivateTest, constantPhiTest); 

};

#endif
