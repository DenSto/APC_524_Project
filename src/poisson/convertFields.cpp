#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "poisson.hpp"

/// derives E from scalar potential phi: E = -grad phi
/*! Makes three calls to phiToESingleComp which performs actual computation
 */ 
void Poisson_Solver::phiToE() { 
    int xdir=0; 
    int ydir=1; 
    int zdir=2; 
    
    phiToESingleComp_(ExID_,xdir); 
    phiToESingleComp_(EyID_,ydir); 
    phiToESingleComp_(EzID_,zdir); 
}; 

/// Calculates a single component of E from phi 
/*! fieldID is a field's fieldID (ExID_, EyID_, or EzID_) \n
 * dir is (0,1,2) for (x,y,z) which must match the component of the fieldID being solved for \n
 */ 
void Poisson_Solver::phiToESingleComp_(const int fieldID, const int dir) { 

    // get the field type and pointer from ID 
    int type = fieldType_[fieldID]; 
    double*** field  = fieldPtr_[fieldID]; 

    // get the last physical point in each direction 
    const int xdir=0; 
    const int ydir=1; 
    const int zdir=2; 
    int iEnd = fieldSize_[type][xdir]; 
    int jEnd = fieldSize_[type][ydir]; 
    int kEnd = fieldSize_[type][zdir]; 

    // get the step size for the derivative
    double h; 
    switch (dir) { 
        case (xdir): h = dx_; break; 
        case (ydir): h = dy_; break; 
        case (zdir): h = dz_; break; 
    } 
    
    // loop over all physical i,j,k to set E = - grad phi
    double nextVal; 
    int i,j,k;  
    for (i=iBeg_; i < iEnd+1; ++i) { 
        for (j=jBeg_; j < jEnd+1; ++j) { 
            for (k=kBeg_; k < kEnd+1; ++k) { 
                switch (dir) { 
                    case (xdir): 
                        nextVal = phi1_[i+1][j][k];
                        break; 
                    case (ydir): 
                        nextVal = phi1_[i][j+1][k];
                        break; 
                    case (zdir): 
                        nextVal = phi1_[i][j][k+1];
                        break; 
                } 
                /* although this exression resembles a forward
                 * difference, it is actually a correct center
                 * difference. E is located a half step between 
                 * the next and prev phi points, and the indexing
                 * of phi and E are also offset */ 
                field[i][j][k] = -(nextVal - phi1_[i][j][k])/h;  
            } 
        } 
    } 
}; 

/// derives A from vector potential A: B = curl A 
/*! Makes three separate calls to AToBSingleComp to perform calculation
 */ 
void Poisson_Solver::AToB() { 
    int xdir=0; 
    int ydir=1; 
    int zdir=2; 
    
    AToBSingleComp_(BxID_,xdir); 
    AToBSingleComp_(ByID_,ydir); 
    AToBSingleComp_(BzID_,zdir); 
}; 

/// calculates a single component of B from A 
/*! fieldID is the field ID of the component to be solved for (BxID_, ByID_, or BzID_) \n 
 * dir is the direction corresonding to the component being solved for
 */ 
void Poisson_Solver::AToBSingleComp_(const int fieldID, const int dir) { 
    
    // get the field type and pointer from ID 
    int type = fieldType_[fieldID]; 
    double*** field  = fieldPtr_[fieldID]; 

    // get the last physical point in each direction 
    const int xdir=0; 
    const int ydir=1; 
    const int zdir=2; 
    int iEnd = fieldSize_[type][xdir]; 
    int jEnd = fieldSize_[type][ydir]; 
    int kEnd = fieldSize_[type][zdir]; 

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
    else if (fieldID == BzID_) { 
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
    for (i=iBeg_; i < iEnd+1; ++i) { 
        for (j=jBeg_; j < jEnd+1; ++j) { 
            for (k=kBeg_; k < kEnd+1; ++k) { 
                switch (dir) { 
                    case (xdir): 
                        nextValk = Ak[i][j+1][k]; 
                        nextValj = Ak[i][j][k+1]; 
                        break; 
                    case (ydir): 
                        nextValk = Ak[i][j][k+1]; 
                        nextValj = Ak[i+1][j][k]; 
                        break; 
                    case (zdir): 
                        nextValk = Ak[i+1][j][k]; 
                        nextValj = Ak[i][j+1][k]; 
                        break; 
                }
                /* although this exression resembles a forward
                 * difference, it is actually a correct center
                 * difference. B is located a half step between 
                 * the next and prev A points, and the indexing
                 * of A and B are also offset */ 
                field[i][j][k] = (nextValk - Ak[i][j][k])/hj + 
                    (nextValj - Aj[i][j][k])/hk;  
            } 
        } 
    } 
}; 
