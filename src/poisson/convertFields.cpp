#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "poisson.hpp"

/// derives E from scalar potential phi: E = -grad phi
void Poisson_Solver::phiToE_() { 
    int xdir=0; 
    int ydir=1; 
    int zdir=2; 
    
    phiToESingleComp_(ExID_,xdir,dx_); 
    phiToESingleComp_(EyID_,ydir,dy_); 
    phiToESingleComp_(EzID_,zdir,dz_); 
}; 


void Poisson_Solver::phiToESingleComp_(const int fieldID, const int dir, const double h) { 

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
                field[i][j][k] = -(nextVal - phi1_[i][j][k])/h;  
            } 
        } 
    } 
}; 

/// derives A from vector potential A: B = curl A 
void Poisson_Solver::AToB_() { 
    printf("Warning: this method currently empty"); 
}; 

