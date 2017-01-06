#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "poisson.hpp"
#include "../globals.hpp"
#include "../domain/domain.hpp"

Poisson_Solver::Poisson_Solver(Domain *domain, Input_Info_t *input_info) :
  Grid(domain->getnxyz(), domain->getnGhosts(), domain->getxyz0(), domain->getLxyz()),
  nFieldsPoisson_(8),
  phi1ID_(13),
  phi2ID_(14),
  Ax1ID_(15),
  Ay1ID_(16),
  Az1ID_(17),
  Ax2ID_(18),
  Ay2ID_(19),
  Az2ID_(20), 
  xdir_(1), 
  ydir_(2),
  zdir_(3)
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

  double sourceMult = -4*M_PI;
  double convTol = .01;
  run_poisson_solver_(phi1ID_,phi1_,phi2_,rho_,convTol,sourceMult);
  phiToE();

  sourceMult = -4*M_PI; 
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

  // directions 
  int xside=1; 
  int yside=2; 
  int zside=3; 
  
  // limits 
  int iEnd = sideToIndex_(xside,fieldID)+1;  
  int jEnd = sideToIndex_(yside,fieldID)+1; 
  int kEnd = sideToIndex_(zside,fieldID)+1; 
  
  //loop Jacobi method until convergence!
  int iternum = -1;
  do {
    iternum++;
    //iterate over entire grid. Note boundary conditions must be supplied!
    for ( int i=iBeg_; i<iEnd; i++ ) {
      for ( int j=jBeg_; j<jEnd; j++ ) {
    	for ( int k=kBeg_; k<kEnd; k++ ) {
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
    domain_->PassFields(this, input_info_,-3,0); // sendID=-3 to bundle phi1(1), phi2(1), A1(3), and A2(3) fieldIDs, op=0 to pass and replace

    //Determine global convergence of jacobi method across all MPI domains
#if USE_MPI
    maxDiff = domain_->GetMaxValueAcrossDomains(maxDiff);
#endif
    if (maxDiff < convergenceTol) jacobi_method_converged = true;
  }while( !jacobi_method_converged );

  //If last calculated pass updated the work array (u1), copy it to the solution array (u0).
  if ( iternum % 2 == 0 ) std::swap(u0,u1);

}

/// returns size of ghost cell data to send
/*! sendID is an integer specifying which fields are intended to be packaged into the ghost vector. -3 for potentials (phi,A), -2 for sources (rho,J), -1 for fields (EB), fieldID for any individual field (e.g. ExID_) \n 
 * It is of length equal to the number of fields being sent times the maximum number of total points in any plane, so that it will be large enough to send the maximum amount of data in a single plane of any of the fields. 
 */ 
int Poisson_Solver::getGhostVecSize(const int sendID) { 
    assert(sendID > -4 && sendID < nFieldsTotal_);
    switch (sendID) {
        case -3: return nFieldsPoisson_*maxPointsInPlane_; break; // A,phi
        case -2: return 4*maxPointsInPlane_; break; // J,rho
        case -1: return 6*maxPointsInPlane_; break; // E,B
        default: return maxPointsInPlane_; break; // single field
    }
}; 

/// bundles the data in the ghost cells to send
/*! side = -/+ 1 for left/right x direction, -/+ 2 for y, -/+ 3 for z \n
 * ghostVec is the vector to store the data in, which must be of length ghostVecSize_ (can be determined with getGhostVecSize) \n
 * sendID = -1 to get EB fields, -2 for rho/J sources, -3 for phi/A potentials, or sendID = an individual field ID (e.g. ExID_) to get just that field (used for Poisson updating for example) \n
 * Stores the data of the E,B,J fields along the specified boundary plane into a 1D array to be sent with a single MPI call. \n 
 * If sendID = -1 (as used in each time step update), stores in order: Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz. \n
 * If sendID = -2 (as used in Poisson iteration), stores in order: phi1,phi2,Ax1,Ay1,Az1,Ax2,Ay2,Az2 \n
 * If sendID = -3 (as used in Poisson initialization), stores in order: Jx,Jy,Jz,rho \n
 * ghostVec can (and should) be unpacked with setGhostVec function 
 */ 
void Poisson_Solver::getGhostVec(const int side, double* ghostVec, int sendID) {
    assert(-4 < sendID && sendID < nFieldsTotal_); 
    
    // create a temporary vector to store slices in 
    int n = maxPointsInPlane_;
    double* tmpVec = sliceTmp_; 

    // offset = 0 to get from the first/last physical cells 
    int offset=0;

    // determine number of fields being sent 
    int nfields; 
    switch (sendID) { 
        case -3: nfields=nFieldsPoisson_; break; 
        case -2: nfields=4; break; 
        case -1: nfields=6; break; 
        default: nfields=1; break; 
    }
    
    // "loop" over all fields to package 
    int begdex; 
    double*** field; 
    int fieldID,ifield;
    for (ifield=0; ifield<nfields; ++ifield) { 
        begdex=ifield*n; 
        switch (sendID) { 
            case -3: // send A/phi 
                switch (ifield) { 
                    case 0: fieldID = phi1ID_; break; 
                    case 1: fieldID = phi2ID_; break; 
                    case 2: fieldID = Ax1ID_; break; 
                    case 3: fieldID = Ay1ID_; break; 
                    case 4: fieldID = Az1ID_; break; 
                    case 5: fieldID = Ax2ID_; break; 
                    case 6: fieldID = Ay2ID_; break; 
                    case 7: fieldID = Az2ID_; break; 
                }; 
                break; 
            case -2: // send J/rho 
                switch (ifield) { 
                    case 0: fieldID = JxID_; break; 
                    case 1: fieldID = JyID_; break; 
                    case 2: fieldID = JzID_; break; 
                    case 3: fieldID = rhoID_; break; 
                }; 
                break; 
            case -1: // send E/B
                switch (ifield) { 
                    case 0: fieldID = ExID_; break; 
                    case 1: fieldID = EyID_; break; 
                    case 2: fieldID = EzID_; break; 
                    case 3: fieldID = BxID_; break; 
                    case 4: fieldID = ByID_; break; 
                    case 5: fieldID = BzID_; break; 
                }
                break; 
            default: fieldID = sendID; break; // send individual field 
        }; 
        field = fieldPtr_[fieldID]; 
        // slice the given field 
        sliceMatToVec_(fieldID,side,offset,tmpVec); 
        // store the slice in ghostVec 
        std::copy(tmpVec,tmpVec + n ,ghostVec + begdex); 
    } 
}; 

/// unbundles the data in the ghost cells to send
/*! side = -/+ 1 for left/right x direction, -/+ 2 for y, -/+ 3 for z \n
 * ghostVec is the vector to read the data from, which must be of length ghostVecSize_ (can be determined with getGhostVecSize) \n
 * sendID = -1 to set JEB fields, or sendID = an individual field ID (e.g. ExID_) to set just that field (used for Poisson updating for example) \n
 * Sets the data of the E,B,J fields along the specified boundary plane from the 1D array ghostVec to be received with a single MPI call. \n 
 * If sendID = -1 (as used in each time step update), fields are read and set in order: Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz. \n
 * If sendID = -2 (as used in Poisson iteration), fields are read and set in order: phi1,phi2,Ax1,Ay1,Az1,Ax2,Ay2,Az2 \n
 * If sendID = -3 (as used in Poisson initialization), stores in order: Jx,Jy,Jz,rho \n
 * ghostVec can (and should) be generated with getGhostVec function 
 */ 
void Poisson_Solver::setGhostVec(const int side, double* ghostVec, int sendID, int op) {
    assert(-4 < sendID && sendID < nFieldsTotal_); 
    
    // create a temporary vector to store slices in 
    int n = maxPointsInPlane_;
    double* tmpVec = sliceTmp_; 
    
    // offset = +1 to set into the RHS ghost vectors
    // offset = -1 to set into the LHS ghost vectors 
    int offset = 1; 
    if (side < 0) { 
        offset = -1; 
    };

    // determine number of fields being sent 
    int nfields; 
    switch (sendID) { 
        case -3: nfields=nFieldsPoisson_; break; 
        case -2: nfields=4; break; 
        case -1: nfields=6; break; 
        default: nfields=1; break; 
    }
    
    // "loop" over all fields to unpackage 
    int begdex,enddex; 
    double*** field; 
    int fieldID,ifield;
    for (ifield=0; ifield<nfields; ++ifield) { 
        begdex=ifield*n; 
        enddex=(ifield+1)*n; 
        // store the relevant portion for ghostVec into tmpVec
        std::copy(ghostVec + begdex, ghostVec + enddex, tmpVec);
        switch (sendID) { 
            case -3: // send A/phi 
                switch (ifield) { 
                    case 0: fieldID = phi1ID_; break; 
                    case 1: fieldID = phi2ID_; break; 
                    case 2: fieldID = Ax1ID_; break; 
                    case 3: fieldID = Ay1ID_; break; 
                    case 4: fieldID = Az1ID_; break; 
                    case 5: fieldID = Ax2ID_; break; 
                    case 6: fieldID = Ay2ID_; break; 
                    case 7: fieldID = Az2ID_; break; 
                }; 
                break; 
            case -2: // send J/rho 
                switch (ifield) { 
                    case 0: fieldID = JxID_; break; 
                    case 1: fieldID = JyID_; break; 
                    case 2: fieldID = JzID_; break; 
                    case 3: fieldID = rhoID_; break; 
                }; 
                break; 
            case -1: // send E/B
                switch (ifield) { 
                    case 0: fieldID = ExID_; break; 
                    case 1: fieldID = EyID_; break; 
                    case 2: fieldID = EzID_; break; 
                    case 3: fieldID = BxID_; break; 
                    case 4: fieldID = ByID_; break; 
                    case 5: fieldID = BzID_; break; 
                }
                break; 
            default: fieldID = sendID; break; // send individual field 
        }; 
        field = fieldPtr_[fieldID]; 
        // unslice the given field 
        unsliceMatToVec_(fieldID,side,offset,tmpVec,op); 
    } 
}; 
