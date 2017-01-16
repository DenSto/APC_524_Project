#include <algorithm>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <assert.h>
#include "poisson.hpp"
#include "../globals.hpp"
//#include "../domain/domain.hpp"

#define CONV_CRIT 1e-9 // absolute convergence criteria
#define CONV_RATIO 1e-6 // relative convergence criteria
#define ITER_MAX 10000 //maximum number of Poisson solver iterations before failure.

Poisson_Solver::Poisson_Solver(Domain *domain, Input_Info_t *input_info) :
  Grid(domain->getnxyz(), 1, domain->getxyz0(), domain->getLxyz()),
  nFieldsPoisson_(8),
  phi1ID_(13),
  Ax1ID_(14),
  Ay1ID_(15),
  Az1ID_(16),
  phi2ID_(17),
  Ax2ID_(18),
  Ay2ID_(19),
  Az2ID_(20), 
  xdir_(1), 
  ydir_(2),
  zdir_(3)
{
  // check field ID is specified in ways that are compatible
  // with the needs of physical boundary condition
  assert(phi2ID_==phi1ID_+4);
  assert(Ax2ID_==Ax1ID_+4);
  assert(Ay2ID_==Ay1ID_+4);
  assert(Az2ID_==Az1ID_+4);

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

  // determine convergence criteria of poissonsolver
  double vmax,vmin;

  std::vector<double> phi;
  for(int i=0;i<6;i++){phi.push_back(input_info_->bound_phi[i]);}
  vmax = *std::max_element(phi.begin(),phi.end());
  vmin = *std::min_element(phi.begin(),phi.end());
  conv_phi_ = vmax-vmin;
  if(debug>2&&rank_MPI==0)fprintf(stderr,"conv_phi_=%f\n",conv_phi_);

  std::vector<double> Ai;
  for(int i=0;i<6;i++){
      Ai.push_back(input_info_->bound_Ax[i]);
      Ai.push_back(input_info_->bound_Ay[i]);
      Ai.push_back(input_info_->bound_Az[i]);
   }
  vmax = *std::max_element(Ai.begin(),Ai.end());
  vmin = *std::min_element(Ai.begin(),Ai.end());
  conv_A_ = vmax-vmin;
  if(debug>2&&rank_MPI==0)fprintf(stderr,"conv_A_=%f\n",conv_A_);

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

//! return fieldID for phi and Ax,Ay,Az
int Poisson_Solver::getFieldID(const std::string &fieldStr){
  int ID;
  if(fieldStr == "phi" || fieldStr == "phi1"){ID = phi1ID_;}
  else if(fieldStr == "Ax" || fieldStr == "Ax1"){ID = Ax1ID_;}
  else if(fieldStr == "Ay" || fieldStr == "Ay1"){ID = Ay1ID_;}
  else if(fieldStr == "Az" || fieldStr == "Az1"){ID = Az1ID_;}
  else if(fieldStr == "phi2"){ID = phi2ID_;}
  else if(fieldStr == "Ax2" ){ID = Ax2ID_;}
  else if(fieldStr == "Ay2" ){ID = Ay2ID_;}
  else if(fieldStr == "Az2" ){ID = Az2ID_;}
  else{ID = Grid::getFieldID(fieldStr);} //call base class function
  return ID;
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
void Poisson_Solver::InitializeFields(Input_Info_t *input_info) {
  //Requires input_info argument to properly inherit virtual Grid::InitializeFields(Input_Info_t *input_info)

  if(rank_MPI==0)printf("        Initializing fields by solving Poisson's equation...\n");

  double sourceMult = -4*M_PI;
  double convTol;

  convTol = std::max(conv_phi_*CONV_RATIO,CONV_CRIT);
  if(debug>1&&rank_MPI==0)printf("        solving for phi...\n"); 
  run_poisson_solver_(phi1ID_,phi2ID_,phi1_,phi2_,rho_,convTol,sourceMult);
  phiToE();

  sourceMult = -4*M_PI; 
  convTol = std::max(conv_A_*CONV_RATIO,CONV_CRIT);
  if(debug>1&&rank_MPI==0)printf("        solving for Ax...\n"); 
  run_poisson_solver_(Ax1ID_,Ax2ID_,Ax1_,Ax2_,Jx_,convTol,sourceMult);

  if(debug>1&&rank_MPI==0)printf("        solving for Ay...\n"); 
  run_poisson_solver_(Ay1ID_,Ay2ID_,Ay1_,Ay2_,Jy_,convTol,sourceMult);
  
  if(debug>1&&rank_MPI==0)printf("        solving for Az...\n"); 
  run_poisson_solver_(Az1ID_,Az2ID_,Az1_,Az2_,Jz_,convTol,sourceMult);
  AToB();

  executeBC(-1); //Execute boundary conditions for E and B (-1) with replace operation (0).
  if(rank_MPI==0)printf("        Poisson initialization complete!\n");
}

void Poisson_Solver::run_poisson_solver_(const int fieldID1, const int fieldID2, double*** u1, double*** u2,double*** R,double convergenceTol,double sourceMult) {
  //fieldID1 is the field u1's ID, e.g. fieldID1 for Ax1_ is Ax1ID_
  //fieldID1 is the field u2's ID, e.g. fieldID2 for Ax2_ is Ax2ID_
  //u1 is the first guess at a solution to Poisson's eqn
  //u2 is the work array of equal size
  //R is the 'source' array for Poisson's eqn
  //convergenceTol is the needed (absolute) solution accuracy

  //Define constants used to iterate Poisson's equation
  double celldist2 = pow(dx_, 2.0)*pow(dy_, 2.0) + pow(dy_, 2.0)*pow(dz_, 2.0) + pow(dx_, 2.0)*pow(dz_, 2.0);
  double ax = pow(dy_, 2.0) * pow(dz_, 2.0) / (2.0 * celldist2);
  double ay = pow(dx_, 2.0) * pow(dz_, 2.0) / (2.0 * celldist2);
  double az = pow(dx_, 2.0) * pow(dy_, 2.0) / (2.0 * celldist2);
  double af = pow(dx_, 2.0) * pow(dy_, 2.0) * pow(dz_, 2.0) / (2.0 * celldist2);

  // limits
  int iEnd = nxReal_ + iBeg_;
  int jEnd = nyReal_ + jBeg_;
  int kEnd = nzReal_ + kBeg_;

  //initialize poisson convergence variables
  bool jacobi_method_converged = false;
  double maxDiff = 0.0;
  double absDiff = 0.0;
  double avgR = 0.0;
  long countR = 0;
  for ( int i=iBeg_; i<iEnd; i++ ) {
    for ( int j=jBeg_; j<jEnd; j++ ) {
      for ( int k=kBeg_; k<kEnd; k++ ) {
	if (R[i][j][k] != 0) {
	  avgR += fabs(R[i][j][k]);
	  countR += 1;
	}
      }
    }
  }
  avgR /= countR;
  convergenceTol = std::min(convergenceTol, fabs(af * avgR * sourceMult ));

  //loop Jacobi method until convergence!
  int iternum = -1;
  do {
    iternum++;

    // supply boundary conditions
    if ( iternum % 2 == 0 ) {
      executeBC(fieldID1);
    } else {
      executeBC(fieldID2);
    }

    //reset poisson convergence variable
    maxDiff = 0.0;

    //iterate over physical grid points only.
    for ( int i=iBeg_; i<iEnd; i++ ) {
      for ( int j=jBeg_; j<jEnd; j++ ) {
	for ( int k=kBeg_; k<kEnd; k++ ) {
	  //Calculate poisson step
	  if ( iternum % 2 == 0 ) {
	    u2[i][j][k] = ax*(u1[i-1][j][k]+u1[i+1][j][k]) 
	      + ay*(u1[i][j-1][k]+u1[i][j+1][k]) 
	      + az*(u1[i][j][k-1]+u1[i][j][k+1]) 
	      - af*R[i][j][k]*sourceMult;
	  } else {
	    u1[i][j][k] = ax*(u2[i-1][j][k]+u2[i+1][j][k]) 
	      + ay*(u2[i][j-1][k]+u2[i][j+1][k]) 
	      + az*(u2[i][j][k-1]+u2[i][j][k+1]) 
	      - af*R[i][j][k]*sourceMult;
	  }
	}
      }
    }

    // calculate error in convergence (after update step above, BC's shouldnt matter)
    for ( int i=iBeg_; i<iEnd; i++ ) {
      for ( int j=jBeg_; j<jEnd; j++ ) {
        for ( int k=kBeg_; k<kEnd; k++ ) {
          //Calculate poisson step
          if ( iternum % 2 == 0 ) {
            absDiff = fabs(u2[i][j][k]-u1[i][j][k]);
          } else {
            absDiff = fabs(u1[i][j][k]-u2[i][j][k]);
          }
          //Retain the largest error
          if (absDiff > maxDiff) {
	    maxDiff = absDiff;
	  }
        }
      }
    }

    //Determine global convergence of jacobi method across all MPI domains
#if USE_MPI
    maxDiff = domain_->GetMaxValueAcrossDomains(maxDiff);
#endif

    if (maxDiff < convergenceTol) jacobi_method_converged = true;

  }while( !jacobi_method_converged && iternum < ITER_MAX);

  if (iternum < ITER_MAX) {
    if (debug) printf("Poisson converged with maxDiff=%e!\n",maxDiff);
  } else {
    printf("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
    printf("WARNING: After %d steps, Poisson FAILED to converge, with maxDiff=%.10e, and convTol=%.10e!\n",ITER_MAX,maxDiff,convergenceTol);
    printf("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
  }

  //If last calculated pass updated the work array (u2), copy it to the solution array (u1).
  if ( iternum % 2 == 0 ) std::swap(u1,u2);

}

/// returns size of ghost cell data to send
/*! sendID is an integer specifying which fields are intended to be packaged into the ghost vector. -3 for potentials (phi,A), -2 for sources (rho,J), -1 for fields (EB), fieldID for any individual field (e.g. ExID_) \n 
 * It is of length equal to the number of fields being sent times the maximum number of total points in any plane, so that it will be large enough to send the maximum amount of data in a single plane of any of the fields. 
 */ 
int Poisson_Solver::getGhostVecSize(const int sendID) { 
    assert(sendID > -4 && sendID < nFieldsTotal_);
    switch (sendID) {
        case -3: return nFieldsPoisson_*maxPointsInPlane_; break; // A,phi
        case -2: return 3*maxPointsInPlane_; break; // J,rho
        case -1: return 6*maxPointsInPlane_; break; // E,B
        default: return maxPointsInPlane_; break; // single field
    }
}; 

/// bundles the data in the ghost cells to send
/*! side = -/+ 1 for left/right x direction, -/+ 2 for y, -/+ 3 for z \n
 * ghostVec is the vector to store the data in, which must be of length ghostVecSize_ \n
 * sendID = -1 to get EB fields, 
 *          -2 for rho/J sources, 
 *          -3 for phi/A potentials, 
 *        = an individual field ID (e.g. ExID_) to get just that field \n
 * ghostVec can (and should) be unpacked with setGhostVec function 
 */ 
void Poisson_Solver::getGhostVec(const int side, double* ghostVec, int sendID) {
    assert(-4 < sendID && sendID < nFieldsTotal_); 
    
    int offset; // offset from physical boundary
    if(sendID==-2){// get Jj,Jk,Rho on right ghost
        assert(side>0);
        offset = 1; // right ghost
    }else{// get physical values
        offset = 0; // physical
    }
    if(debug>1)fprintf(stderr,"rank=%d:get on side=%d with offset=%d,for field=%d\n",
                                  rank_MPI,side,offset,sendID);

    // create a temporary vector to store slices in 
    int n = maxPointsInPlane_;
    double* tmpVec = sliceTmp_; 

    const int xside = 1; 
    // const int yside = 2; // not used 
    const int zside = 3; 
    
    // determine number of fields being sent 
    int nfields; 
    switch (sendID) { 
        case -3: nfields=nFieldsPoisson_; break; 
        case -2: nfields=3; break; 
        case -1: nfields=6; break; 
        default: nfields=1; break; 
    }
    
    // "loop" over all fields to package 
    int begdex; 
    double*** field; 
    int fieldID; 
    int ifield;
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
                    // determine the 2 orthogonal J components to send 
                    case 0: 
                        if (abs(side) == xside) { 
                            fieldID = JyID_; 
                        } else { 
                            fieldID = JxID_; 
                        } 
                        break; 
                    case 1: 
                        if (abs(side) == zside) { 
                            fieldID = JyID_; 
                        } else { 
                            fieldID = JzID_; 
                        } 
                        break; 
                    case 2: fieldID = rhoID_; break; 
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
        
        // slice the given field with appropriate offset 
        sliceMatToVec_(fieldID,side,offset,tmpVec); 
        // store the slice in ghostVec 
        std::copy(tmpVec,tmpVec + n ,ghostVec + begdex); 
    } 
}; 

/// unbundles the data in the ghost cells to send
/*! side = -/+ i for left/right i-th direction \n
 * ghostVec is the vector to read the data from, which must be of length ghostVecSize_ \n
 * sendID = -3 phi and A fields
 *          -2 to set Jj,Jk, rho fields
 *          -1 to set EB fields 
 *          an individual field ID (e.g. ExID_) to set just that field \n
 * ghostVec can (and should) be generated with getGhostVec function 
 */
void Poisson_Solver::setGhostVec(const int side, double* ghostVec, int sendID) {
    assert(-4 < sendID && sendID < nFieldsTotal_); 
    
    int offset; // offset from physical boundary
    int option; // for unslice,0:replace, 1:sum
    if(sendID==-2){// sum Rho and Jj,Jk to left physical
        assert(side<0); //left
        offset = 0; // physical
        option = 1; // sum
    }else{// set ghost values by replace
        offset = abs(side)/side; // ghost
        option = 0; // replace
    }
    if(debug>1)fprintf(stderr,"rank=%d:set side=%d, offset=%d, field=%d, option=%d\n",
                                  rank_MPI,side,offset,sendID,option);

    // create a temporary vector to store slices in 
    int n = maxPointsInPlane_;
    double* tmpVec = sliceTmp_; 
    
    const int xside = 1; 
    // const int yside = 2; // not used 
    const int zside = 3; 
    
    // determine number of fields being sent 
    int nfields; 
    switch (sendID) { 
        case -3: nfields=nFieldsPoisson_; break; 
        case -2: nfields=3; break; 
        case -1: nfields=6; break; 
        default: nfields=1; break; 
    }
    
    // "loop" over all fields to unpackage 
    int begdex,enddex; 
    double*** field; 
    int fieldID; 
    int ifield;
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
            case -2: // set J/rho 
                switch (ifield) { 
                    // determine the 2 orthogonal J components to set 
                    case 0: 
                        if (abs(side) == xside) { 
                            fieldID = JyID_; 
                        } else { 
                            fieldID = JxID_; 
                        } 
                        break; 
                    case 1: 
                        if (abs(side) == zside) { 
                            fieldID = JyID_; 
                        } else { 
                            fieldID = JzID_; 
                        } 
                        break; 
                    case 2: fieldID = rhoID_; break; 
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
        unsliceMatToVec_(fieldID,side,offset,tmpVec,option); 
    } 
}; 
