#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "poisson.hpp"

/// derives E from scalar potential phi: E = -grad phi
/*! Makes three calls to phiToESingleComp which performs actual computation
 */ 
void Poisson_Solver::phiToE() { 
    phiToESingleComp_(ExID_,xdir_); 
    phiToESingleComp_(EyID_,ydir_); 
    phiToESingleComp_(EzID_,zdir_); 
}; 

/// Calculates a single component of E from phi 
/*! fieldID is a field's fieldID (ExID_, EyID_, or EzID_) \n
 * dir is (0,1,2) for (x,y,z) which must match the component of the fieldID being solved for \n
 */ 
void Poisson_Solver::phiToESingleComp_(const int fieldID, const int dir) { 
    assert (0 < dir && dir < 4); 
    assert(fieldID==ExID_ || fieldID==EyID_ || fieldID==EzID_);
    assert(dir==xdir_ || dir==ydir_ || dir==zdir_);

    // get the field pointer from ID 
    double*** field  = fieldPtr_[fieldID]; 

    // limits 
    int iEnd = sideToIndex_(xdir_,fieldID)+1;  
    int jEnd = sideToIndex_(ydir_,fieldID)+1; 
    int kEnd = sideToIndex_(zdir_,fieldID)+1; 

    // get the step size for the derivative
    double h;
    if (dir == xdir_) { 
        h = dx_; 
    } else if (dir == ydir_) { 
        h = dy_; 
    //} else if (dir == zdir_) { 
    } else {
        h = dz_; 
    }; 

    // loop over all physical i,j,k to set E = - grad phi
    double nextVal; 
    int i,j,k;
    for (i=iBeg_; i < iEnd; ++i) {
        for (j=jBeg_; j < jEnd; ++j) {
            for (k=kBeg_; k < kEnd; ++k) {
                if (dir == xdir_) { 
                    nextVal = phi1_[i+1][j][k]; 
                } else if (dir == ydir_) { 
                    nextVal = phi1_[i][j+1][k]; 
                //} else if (dir == zdir_) { 
                } else {
                    nextVal = phi1_[i][j][k+1]; 
                }
                /* although this expression resembles a forward
                 * difference, it is actually a correct center
                 * difference. E is located a half step between 
                 * the next and prev phi points, and the indexing
                 * of phi and E are also offset */
                // - sign crucial for E = - grad phi
                field[i][j][k] = -(nextVal - phi1_[i][j][k])/h;
            }
        }
    }
}; 

/// derives A from vector potential A: B = curl A 
/*! Makes three separate calls to AToBSingleComp to perform calculation
 */ 
void Poisson_Solver::AToB() { 
    AToBSingleComp_(BxID_,xdir_); 
    AToBSingleComp_(ByID_,ydir_); 
    AToBSingleComp_(BzID_,zdir_); 
}; 

/// calculates a single component of B from A 
/*! fieldID is the field ID of the component to be solved for (BxID_, ByID_, or BzID_) \n 
 * dir is the direction corresonding to the component being solved for
 */ 
void Poisson_Solver::AToBSingleComp_(const int fieldID, const int dir) { 
    assert(0 < dir && dir < 4); 
    assert(fieldID==BxID_ || fieldID==ByID_ || fieldID==BzID_);
    assert(dir==xdir_ || dir==ydir_ || dir==zdir_);

    // get the field pointer from ID 
    double*** field  = fieldPtr_[fieldID]; 
    
    // limits 
    int iEnd = sideToIndex_(xdir_,fieldID)+1;  
    int jEnd = sideToIndex_(ydir_,fieldID)+1; 
    int kEnd = sideToIndex_(zdir_,fieldID)+1; 

    // set the two components of A to use in curl 
    double*** Ak;
    double*** Aj; 
    double hk,hj; 
    if (fieldID == BxID_) { 
        Ak = Az1_; 
        hk = dz_; 
        Aj = Ay1_; 
        hj = dy_; 
    } 
    else if (fieldID == ByID_) { 
        Ak = Ax1_; 
        hk = dx_; 
        Aj = Az1_; 
        hj = dz_; 
    } 
    //else if (fieldID == BzID_) { 
    else{
        Ak = Ay1_; 
        hk = dy_; 
        Aj = Ax1_; 
        hj = dx_; 
    };

    // loop over all physical i,j,k to set B = curl A
    // indices are named j,k imagining that the component to be 
    // solved or is B_i so that 
    // B_i = epsilon_ijk partial_j A_k, or expanded
    // B_i = d_j A_k - d_k A_j, where 
    // i x j = k is the coordinate system
    double nextValk,nextValj; 
    int i,j,k;  
    for (i=iBeg_; i < iEnd; ++i) { 
        for (j=jBeg_; j < jEnd; ++j) { 
            for (k=kBeg_; k < kEnd; ++k) { 
                if (dir == xdir_) { 
                    nextValk = Ak[i][j+1][k]; 
                    nextValj = Ak[i][j][k+1]; 
                } else if (dir == ydir_) { 
                    nextValk = Ak[i][j][k+1]; 
                    nextValj = Ak[i+1][j][k]; 
                //} else if (dir == zdir_) { 
                } else {
                    nextValk = Ak[i+1][j][k]; 
                    nextValj = Ak[i][j+1][k]; 
                }
                /* although this expression resembles a forward
                 * difference, it is actually a correct center
                 * difference. B is located a half step between 
                 * the next and prev A points, and the indexing
                 * of A and B are also offset */ 
                // hj and hk intentionally do not match Ak and Aj
                // that they are dividing
                field[i][j][k] = (nextValk - Ak[i][j][k])/hj + 
                    (nextValj - Aj[i][j][k])/hk;  
            } 
        } 
    } 
}; 
