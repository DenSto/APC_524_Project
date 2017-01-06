#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "poisson.hpp"
#include "../globals.hpp"
//#include "../domain/domain.hpp"

Poisson_Solver::Poisson_Solver(Domain *domain, Input_Info_t *input_info) :
  Grid(domain->getnxyz(), domain->getnGhosts(), domain->getxyz0(), domain->getLxyz()),
  phi1ID_(13),
  phi2ID_(14),
  Ax1ID_(15),
  Ay1ID_(16),
  Az1ID_(17),
  Ax2ID_(18),
  Ay2ID_(19),
  Az2ID_(20)
{
  phi1_=newField_(phi1ID_);
  phi2_=newField_(phi2ID_);
  Ax1_=newField_(Ax1ID_);
  Ay1_=newField_(Ay1ID_);
  Az1_=newField_(Az1ID_);
  Ax2_=newField_(Ax2ID_);
  Ay2_=newField_(Ay2ID_);
  Az2_=newField_(Az2ID_);

  setPoissonFieldType_();
  setPoissonFieldPtr_();

  domain_ = domain;
  input_info_ = input_info;
}

Poisson_Solver::~Poisson_Solver() {
  deleteField_(phi1_,phi1ID_);
  deleteField_(phi2_,phi2ID_);
  deleteField_(Ax1_,Ax1ID_);
  deleteField_(Ay1_,Ay1ID_);
  deleteField_(Az1_,Az1ID_);
  deleteField_(Ax2_,Ax2ID_);
  deleteField_(Ay2_,Ay2ID_);
  deleteField_(Az2_,Az2ID_);
}

/// Same as Grid::setFieldType_ for phi,A arrays unique to Poisson
void Poisson_Solver::setPoissonFieldType_() { 
    fieldType_[phi1ID_]=vertID_; 
    fieldType_[phi2ID_]=vertID_; 
    fieldType_[Ax1ID_]=edgeXID_; 
    fieldType_[Ay1ID_]=edgeYID_; 
    fieldType_[Az1ID_]=edgeZID_; 
    fieldType_[Ax2ID_]=edgeXID_; 
    fieldType_[Ay2ID_]=edgeYID_; 
    fieldType_[Az2ID_]=edgeZID_; 
}; 

/// Same as Grid::setFieldPtr_ for phi,A arrays unique to Poisson
void Poisson_Solver::setPoissonFieldPtr_() { 
    fieldPtr_[phi1ID_]=phi1_; 
    fieldPtr_[phi2ID_]=phi2_; 
    fieldPtr_[Ax1ID_]=Ax1_; 
    fieldPtr_[Ay1ID_]=Ay1_; 
    fieldPtr_[Az1ID_]=Az1_; 
    fieldPtr_[Ax2ID_]=Ax2_; 
    fieldPtr_[Ay2ID_]=Ay2_; 
    fieldPtr_[Az2ID_]=Az2_; 
};

/// Set vector potential to constant values
void Poisson_Solver::constA(const double vx, const double vy, const double vz) { 
    constField_(Ax1ID_,vx); 
    constField_(Ay1ID_,vy); 
    constField_(Az1ID_,vz); 
    constField_(Ax2ID_,vx); 
    constField_(Ay2ID_,vy); 
    constField_(Az2ID_,vz); 
} 

/// Set scalar potential to constant value
void Poisson_Solver::constPhi(const double v) { 
    constField_(phi1ID_,v); 
    constField_(phi2ID_,v); 
} 

//void Poisson_Solver::initialize_poisson_fields() {
void Poisson_Solver::InitializeFields() {

  if(rank_MPI==0)printf("        Initializing fields by solving Poisson's equation...\n");

  double sourceMult = -4*3.1415926535898;
  double convTol = .01;
  run_poisson_solver_(phi1ID_,phi1_,phi2_,rho_,convTol,sourceMult);
  phiToE();

  sourceMult = -4*3.1415926535898;
  convTol = .1;
  run_poisson_solver_(Ax1ID_,Ax1_,Ax2_,Jx_,convTol,sourceMult);
  run_poisson_solver_(Ay1ID_,Ay1_,Ay2_,Jy_,convTol,sourceMult);
  run_poisson_solver_(Az1ID_,Az1_,Az2_,Jz_,convTol,sourceMult);
  AToB();
}

void Poisson_Solver::run_poisson_solver_(const int fieldID, double*** u0, double*** u1,double*** R,double convergenceTol,double sourceMult) {
  //fieldID is the field's ID, e.g. fieldID for Ax_ is Ax1ID_
  //u0 is the first guess at a solution to Poisson's eqn
  //u1 is the work array of equal size
  //R is the 'source' array for Poisson's eqn
  //convergenceTol is the needed (absolute) solution accuracy

  //Define constants used to iterate Poisson's equation
  double celldist2 = pow(dx_, 2.0) + pow(dy_, 2.0) + pow(dz_, 2.0);
  double ax = pow(dy_, 2.0) * pow(dz_, 2.0) / (2.0 * celldist2);
  double ay = pow(dx_, 2.0) * pow(dz_, 2.0) / (2.0 * celldist2);
  double az = pow(dx_, 2.0) * pow(dy_, 2.0) / (2.0 * celldist2);
  double af = pow(dx_, 2.0) * pow(dy_, 2.0) * pow(dz_, 2.0) / (2.0 * celldist2);

  //initialize iteration variables
  bool jacobi_method_converged = false;
  double maxDiff = 0.0;
  double absDiff = 0.0;

  //loop Jacobi method until convergence!
  int iternum = -1;
  do {
    iternum++;
    //iterate over entire grid. Note boundary conditions must be supplied!
    for ( int i=iBeg_; i<nx_; i++ ) {
      for ( int j=jBeg_; j<ny_; j++ ) {
    	for ( int k=kBeg_; k<nz_; k++ ) {
          if ( iternum % 2 == 0 ) {
            u1[i][j][k] = ax*(u0[i-1][j][k]+u0[i+1][j][k]) + ay*(u0[i][j-1][k]+u0[i][j+1][k]) + 
              az*(u0[i][j][k-1]+u0[i][j][k+1]) - af*R[i][j][k]*sourceMult;
          } else {
            u0[i][j][k] = ax*(u1[i-1][j][k]+u1[i+1][j][k]) + ay*(u1[i][j-1][k]+u1[i][j+1][k]) +
                  az*(u1[i][j][k-1]+u1[i][j][k+1]) - af*R[i][j][k]*sourceMult;
          }

          absDiff = fabs(u1[i][j][k] - u0[i][j][k]);
          if (absDiff > maxDiff) maxDiff = absDiff;
        }
      }
    }

    //Pass fields' data via MPI.
    domain_->PassFields(this, input_info_, phi1ID_);
    domain_->PassFields(this, input_info_, phi2ID_);
    domain_->PassFields(this, input_info_, Ax1ID_);
    domain_->PassFields(this, input_info_, Ay1ID_);
    domain_->PassFields(this, input_info_, Az1ID_);
    domain_->PassFields(this, input_info_, Ax2ID_);
    domain_->PassFields(this, input_info_, Ay2ID_);
    domain_->PassFields(this, input_info_, Az2ID_);

    //Determine global convergence of jacobi method across all MPI domains
#if USE_MPI
    maxDiff = domain_->GetMaxValueAcrossDomains(maxDiff);
#endif
    if (maxDiff < convergenceTol) jacobi_method_converged = true;
  }while( !jacobi_method_converged );

  //If last calculated pass updated the work array (u1), copy it to the solution array (u0).
  if ( iternum % 2 == 0 ) std::swap(u0,u1);

}
