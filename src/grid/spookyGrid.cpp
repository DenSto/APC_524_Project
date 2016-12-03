#include <stdlib.h>
#include <assert.h>
#include <math.h> 
#include "grid.hpp"

// slices a physical plane in the x direction (excludes ghosts) 
void Grid::sliceMatToVec_(double *** const mat, const int side) { 
    int i = sideToIndex_(side); 
    int j,k; // iterators
    int iter=-1;
    for (j=jBeg_; j<jEnd_; ++j) { 
        for (k=kBeg_; k<kEnd_; ++k) { 
            // equivalently: iter = j*nz_ + k
            sliceTmp_[++iter] = mat[i][j][k]; 
        } 
    } 
};

// unslices a physical plane in the x direction (excludes ghosts) 
void Grid::unsliceMatToVec_(double*** mat, const int side) { 
    int i = sideToIndex_(side); 
    int j,k; // iterators
    int iter = -1; 
    for (j=jBeg_; j<jEnd_; ++j) { 
        for (k=kBeg_; k<kEnd_; ++k) { 
            // equivalently: iter = j*nz_ + k
            mat[i][j][k] = sliceTmp_[++iter]; 
        } 
    } 
}; 

// updates ghost cells after evolving the field on physical points
// does this take a performance hit due to cache hit/miss?
void Grid::updateGhostCells() { 
    int i,j,k; // iterators 
    
    // update ghost cells in x direction 
    // iterates over yz plane
    for (j=jBeg_; j<jEnd_; ++j) { 
        for (k=kBeg_; k<kEnd_; ++k) { 
            Ex_[iBeg_-1][j][k]=Ex_[iBeg_][j][k]; 
            Ex_[iEnd_+1][j][k]=Ex_[iEnd_][j][k]; 
            Ey_[iBeg_-1][j][k]=Ey_[iBeg_][j][k]; 
            Ey_[iEnd_+1][j][k]=Ey_[iEnd_][j][k]; 
            Ez_[iBeg_-1][j][k]=Ez_[iBeg_][j][k]; 
            Ez_[iEnd_+1][j][k]=Ez_[iEnd_][j][k]; 

            Bx_[iBeg_-1][j][k]=Bx_[iBeg_][j][k]; 
            Bx_[iEnd_+1][j][k]=Bx_[iEnd_][j][k]; 
            By_[iBeg_-1][j][k]=By_[iBeg_][j][k]; 
            By_[iEnd_+1][j][k]=By_[iEnd_][j][k]; 
            Bz_[iBeg_-1][j][k]=Bz_[iBeg_][j][k]; 
            Bz_[iEnd_+1][j][k]=Bz_[iEnd_][j][k]; 

            Jx_[iBeg_-1][j][k]=Jx_[iBeg_][j][k]; 
            Jx_[iEnd_+1][j][k]=Jx_[iEnd_][j][k]; 
            Jy_[iBeg_-1][j][k]=Jy_[iBeg_][j][k]; 
            Jy_[iEnd_+1][j][k]=Jy_[iEnd_][j][k]; 
            Jz_[iBeg_-1][j][k]=Jz_[iBeg_][j][k]; 
            Jz_[iEnd_+1][j][k]=Jz_[iEnd_][j][k]; 
        } 
    }

    // updates ghost cells in y direction 
    // iterates over xz plane
    for (i=iBeg_; i<iEnd_; ++i) { 
        for (k=kBeg_; k<kEnd_; ++k) { 
            Ex_[i][jBeg_-1][k]=Ex_[i][jBeg_][k]; 
            Ex_[i][jEnd_+1][k]=Ex_[i][jEnd_][k]; 
            Ey_[i][jBeg_-1][k]=Ey_[i][jBeg_][k]; 
            Ey_[i][jEnd_+1][k]=Ey_[i][jEnd_][k]; 
            Ez_[i][jBeg_-1][k]=Ez_[i][jBeg_][k]; 
            Ez_[i][jEnd_+1][k]=Ez_[i][jEnd_][k]; 

            Bx_[i][jBeg_-1][k]=Bx_[i][jBeg_][k]; 
            Bx_[i][jEnd_+1][k]=Bx_[i][jEnd_][k]; 
            By_[i][jBeg_-1][k]=By_[i][jBeg_][k]; 
            By_[i][jEnd_+1][k]=By_[i][jEnd_][k]; 
            Bz_[i][jBeg_-1][k]=Bz_[i][jBeg_][k]; 
            Bz_[i][jEnd_+1][k]=Bz_[i][jEnd_][k]; 

            Jx_[i][jBeg_-1][k]=Jx_[i][jBeg_][k]; 
            Jx_[i][jEnd_+1][k]=Jx_[i][jEnd_][k]; 
            Jy_[i][jBeg_-1][k]=Jy_[i][jBeg_][k]; 
            Jy_[i][jEnd_+1][k]=Jy_[i][jEnd_][k]; 
            Jz_[i][jBeg_-1][k]=Jz_[i][jBeg_][k]; 
            Jz_[i][jEnd_+1][k]=Jz_[i][jEnd_][k]; 
        } 
    }

    // updates ghost cells in z direction 
    // iterates over xy plane 
    for (i=iBeg_; i<iEnd_; ++i) { 
        for (j=jBeg_; j<jEnd_; ++j) { 
            Ex_[i][j][kBeg_-1]=Ex_[i][j][kBeg_]; 
            Ex_[i][j][kEnd_+1]=Ex_[i][j][kEnd_]; 
            Ey_[i][j][kBeg_-1]=Ey_[i][j][kBeg_]; 
            Ey_[i][j][kEnd_+1]=Ey_[i][j][kEnd_]; 
            Ez_[i][j][kBeg_-1]=Ez_[i][j][kBeg_]; 
            Ez_[i][j][kEnd_+1]=Ez_[i][j][kEnd_]; 

            Bx_[i][j][kBeg_-1]=Bx_[i][j][kBeg_]; 
            Bx_[i][j][kEnd_+1]=Bx_[i][j][kEnd_]; 
            By_[i][j][kBeg_-1]=By_[i][j][kBeg_]; 
            By_[i][j][kEnd_+1]=By_[i][j][kEnd_]; 
            Bz_[i][j][kBeg_-1]=Bz_[i][j][kBeg_]; 
            Bz_[i][j][kEnd_+1]=Bz_[i][j][kEnd_]; 

            Jx_[i][j][kBeg_-1]=Jx_[i][j][kBeg_]; 
            Jx_[i][j][kEnd_+1]=Jx_[i][j][kEnd_]; 
            Jy_[i][j][kBeg_-1]=Jy_[i][j][kBeg_]; 
            Jy_[i][j][kEnd_+1]=Jy_[i][j][kEnd_]; 
            Jz_[i][j][kBeg_-1]=Jz_[i][j][kBeg_]; 
            Jz_[i][j][kEnd_+1]=Jz_[i][j][kEnd_]; 
        } 
    }

};

// returns size of ghost cell data to send 
int Grid::getGhostVecSize() { 
    // must bundle 9 scalar fields
    // each of size of the physical zy plane
    return ghostVecSize_; 
}; 

// bundles the data in the ghost cells to send
// sends (in order): Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz
// side is -1 for left side of cell, +1 for right 
// ghostVec is a 1D array storing doubles 
// 9 loops, one for each field (accesses each field sequentially
// better performance due to cache bonuses?
void Grid::getGhostVec(const int side, double* ghostVec) {
    int ncp = nx_-2*nGhosts_; 
    // this is the price of only 3 dimensions on field storage
    sliceMatToVec_(Ex_,side); 
    std::copy(sliceTmp_,sliceTmp_+ncp,ghostVec+(0 * ncp)); 

    sliceMatToVec_(Ey_,side); 
    std::copy(sliceTmp_,sliceTmp_+ncp,ghostVec+(1 * ncp)); 

    sliceMatToVec_(Ez_,side); 
    std::copy(sliceTmp_,sliceTmp_+ncp,ghostVec+(2 * ncp)); 

    sliceMatToVec_(Bx_,side); 
    std::copy(sliceTmp_,sliceTmp_+ncp,ghostVec+(3 * ncp)); 

    sliceMatToVec_(By_,side); 
    std::copy(sliceTmp_,sliceTmp_+ncp,ghostVec+(4 * ncp)); 

    sliceMatToVec_(Bz_,side); 
    std::copy(sliceTmp_,sliceTmp_+ncp,ghostVec+(5 * ncp)); 

    sliceMatToVec_(Jx_,side); 
    std::copy(sliceTmp_,sliceTmp_+ncp,ghostVec+(6 * ncp)); 

    sliceMatToVec_(Jy_,side); 
    std::copy(sliceTmp_,sliceTmp_+ncp,ghostVec+(7 * ncp)); 

    sliceMatToVec_(Jz_,side); 
    std::copy(sliceTmp_,sliceTmp_+ncp,ghostVec+(8 * ncp)); 
}; 

// alternate implementation of getGhostVec
// single loop, pulls from all fields at once
// worse performance due to more cache misses or negligible? 
void Grid::getGhostVecAlt(const int side, double* ghostVec) {
    int realPts = (ny_ - 2*nGhosts_)*(nz_ - 2*nGhosts_); 
    int i = sideToIndex_(side); 
    int j,k; // iterators
    int iter = -1; 
    for (j=jBeg_; j<jEnd_; ++j) { 
        for (k=kBeg_; k<kEnd_; ++k) { 
            ++iter; 
            ghostVec[0*realPts + iter]=Ex_[i][j][k];
            ghostVec[1*realPts + iter]=Ey_[i][j][k];
            ghostVec[2*realPts + iter]=Ez_[i][j][k];
            ghostVec[3*realPts + iter]=Bx_[i][j][k];
            ghostVec[4*realPts + iter]=By_[i][j][k];
            ghostVec[5*realPts + iter]=Bz_[i][j][k];
            ghostVec[6*realPts + iter]=Jx_[i][j][k];
            ghostVec[7*realPts + iter]=Jy_[i][j][k];
            ghostVec[8*realPts + iter]=Jz_[i][j][k];
        } 
    } 
};

// unbundles the data and puts it in the field 
// one loop for each field, aims to minimize cache misses
void Grid::setGhostVec(const int side, const double* ghostVec) { 
    int ncp = nx_-2*nGhosts_; 
    // this is the price of only 3 dimensions on field storage
    std::copy(ghostVec+(0 * ncp),ghostVec+1*ncp -1,sliceTmp_); 
    unsliceMatToVec_(Ex_,side); 

    std::copy(ghostVec+(1 * ncp),ghostVec+2*ncp -1,sliceTmp_); 
    unsliceMatToVec_(Ex_,side); 

    std::copy(ghostVec+(2 * ncp),ghostVec+3*ncp -1,sliceTmp_); 
    unsliceMatToVec_(Ex_,side); 

    std::copy(ghostVec+(3 * ncp),ghostVec+4*ncp -1,sliceTmp_); 
    unsliceMatToVec_(Ex_,side); 

    std::copy(ghostVec+(4 * ncp),ghostVec+5*ncp -1,sliceTmp_); 
    unsliceMatToVec_(Ex_,side); 

    std::copy(ghostVec+(5 * ncp),ghostVec+6*ncp -1,sliceTmp_); 
    unsliceMatToVec_(Ex_,side); 

    std::copy(ghostVec+(6 * ncp),ghostVec+7*ncp -1,sliceTmp_); 
    unsliceMatToVec_(Ex_,side); 

    std::copy(ghostVec+(7 * ncp),ghostVec+8*ncp -1,sliceTmp_); 
    unsliceMatToVec_(Ex_,side); 

    std::copy(ghostVec+(8 * ncp),ghostVec+9*ncp -1,sliceTmp_); 
    unsliceMatToVec_(Ex_,side); 

};

// alternate implementation of setGhostVec
// single loop, sets all fields at once
// worse performance due to more cache misses or negligible? 
void Grid::setGhostVecAlt(const int side, const double* ghostVec) {
    int realPts = (ny_ - 2*nGhosts_)*(nz_ - 2*nGhosts_); 
    int i = sideToIndex_(side);
    int j,k; // iterators
    int iter = -1; 
    for (j=jBeg_; j<jEnd_; ++j) { 
        for (k=kBeg_; k<kEnd_; ++k) { 
            ++iter; 
            Ex_[i][j][k]=ghostVec[0*realPts + iter];
            Ey_[i][j][k]=ghostVec[1*realPts + iter];
            Ez_[i][j][k]=ghostVec[2*realPts + iter];
            Bx_[i][j][k]=ghostVec[3*realPts + iter];
            By_[i][j][k]=ghostVec[4*realPts + iter];
            Bz_[i][j][k]=ghostVec[5*realPts + iter];
            Jx_[i][j][k]=ghostVec[6*realPts + iter];
            Jy_[i][j][k]=ghostVec[7*realPts + iter];
            Jz_[i][j][k]=ghostVec[8*realPts + iter];
        } 
    } 
};

// function to convert -/+ 1 left/right side indicator to index in x direction 
int Grid::sideToIndex_(const int side) { 
    int i; 
    if (side < 0) { 
        i = 1; 
    } 
    else { 
        i = nx_ - (1 + nGhosts_); 
    } 
    return i; 
}; 
