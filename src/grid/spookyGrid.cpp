#include <stdlib.h>
#include <iterator> // used for std::copy
#include <assert.h> // plans to use?
#include <math.h> // unneeded?
#include "grid.hpp"

// slices a physical plane in the x direction (excludes ghosts) 
void Grid::sliceMatToVec_(double *** const mat, const int side) { 
    int i = sideToIndex_(side); 
    int j,k; // iterators
    int iter=-1;
    for (j=jBeg_; j<jEnd_+1; ++j) { 
        for (k=kBeg_; k<kEnd_+1; ++k) { 
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
    for (j=jBeg_; j<jEnd_+1; ++j) { 
        for (k=kBeg_; k<kEnd_+1; ++k) { 
            // equivalently: iter = j*nz_ + k
            mat[i][j][k] = sliceTmp_[++iter]; 
        } 
    } 
}; 

// updates ghost cells after evolving the field on physical points
// does this take a performance hit due to cache hit/miss?
// this method currently requires nGhosts_=1
void Grid::updateGhostCells() { 
    int i,j,k; // iterators 
    int iGhostLeft=iBeg_-1; 
    int iGhostRight=iEnd_+1; 
    int jGhostLeft=jBeg_-1; 
    int jGhostRight=jEnd_+1; 
    int kGhostLeft=kBeg_-1; 
    int kGhostRight=kEnd_+1; 
    
    // update ghost cells in x direction 
    // iterates over yz plane
    for (j=jBeg_; j<jEnd_+1; ++j) { 
        for (k=kBeg_; k<kEnd_+1; ++k) { 
            Ex_[iGhostLeft][j][k]=Ex_[iBeg_][j][k]; 
            Ex_[iGhostRight][j][k]=Ex_[iEnd_][j][k]; 
            Ey_[iGhostLeft][j][k]=Ey_[iBeg_][j][k]; 
            Ey_[iGhostRight][j][k]=Ey_[iEnd_][j][k]; 
            Ez_[iGhostLeft][j][k]=Ez_[iBeg_][j][k]; 
            Ez_[iGhostRight][j][k]=Ez_[iEnd_][j][k]; 

            Bx_[iGhostLeft][j][k]=Bx_[iBeg_][j][k]; 
            Bx_[iGhostRight][j][k]=Bx_[iEnd_][j][k]; 
            By_[iGhostLeft][j][k]=By_[iBeg_][j][k]; 
            By_[iGhostRight][j][k]=By_[iEnd_][j][k]; 
            Bz_[iGhostLeft][j][k]=Bz_[iBeg_][j][k]; 
            Bz_[iGhostRight][j][k]=Bz_[iEnd_][j][k]; 

            Jx_[iGhostLeft][j][k]=Jx_[iBeg_][j][k]; 
            Jx_[iGhostRight][j][k]=Jx_[iEnd_][j][k]; 
            Jy_[iGhostLeft][j][k]=Jy_[iBeg_][j][k]; 
            Jy_[iGhostRight][j][k]=Jy_[iEnd_][j][k]; 
            Jz_[iGhostLeft][j][k]=Jz_[iBeg_][j][k]; 
            Jz_[iGhostRight][j][k]=Jz_[iEnd_][j][k]; 
        } 
    }

    // updates ghost cells in y direction 
    // iterates over xz plane
    for (i=iBeg_; i<iEnd_+1; ++i) { 
        for (k=kBeg_; k<kEnd_+1; ++k) { 
            Ex_[i][jGhostLeft][k]=Ex_[i][jBeg_][k]; 
            Ex_[i][jGhostRight][k]=Ex_[i][jEnd_][k]; 
            Ey_[i][jGhostLeft][k]=Ey_[i][jBeg_][k]; 
            Ey_[i][jGhostRight][k]=Ey_[i][jEnd_][k]; 
            Ez_[i][jGhostLeft][k]=Ez_[i][jBeg_][k]; 
            Ez_[i][jGhostRight][k]=Ez_[i][jEnd_][k]; 

            Bx_[i][jGhostLeft][k]=Bx_[i][jBeg_][k]; 
            Bx_[i][jGhostRight][k]=Bx_[i][jEnd_][k]; 
            By_[i][jGhostLeft][k]=By_[i][jBeg_][k]; 
            By_[i][jGhostRight][k]=By_[i][jEnd_][k]; 
            Bz_[i][jGhostLeft][k]=Bz_[i][jBeg_][k]; 
            Bz_[i][jGhostRight][k]=Bz_[i][jEnd_][k]; 

            Jx_[i][jGhostLeft][k]=Jx_[i][jBeg_][k]; 
            Jx_[i][jGhostRight][k]=Jx_[i][jEnd_][k]; 
            Jy_[i][jGhostLeft][k]=Jy_[i][jBeg_][k]; 
            Jy_[i][jGhostRight][k]=Jy_[i][jEnd_][k]; 
            Jz_[i][jGhostLeft][k]=Jz_[i][jBeg_][k]; 
            Jz_[i][jGhostRight][k]=Jz_[i][jEnd_][k]; 
        } 
    }

    // updates ghost cells in z direction 
    // iterates over xy plane 
    for (i=iBeg_; i<iEnd_+1; ++i) { 
        for (j=jBeg_; j<jEnd_+1; ++j) { 
            Ex_[i][j][kGhostLeft]=Ex_[i][j][kBeg_]; 
            Ex_[i][j][kGhostRight]=Ex_[i][j][kEnd_]; 
            Ey_[i][j][kGhostLeft]=Ey_[i][j][kBeg_]; 
            Ey_[i][j][kGhostRight]=Ey_[i][j][kEnd_]; 
            Ez_[i][j][kGhostLeft]=Ez_[i][j][kBeg_]; 
            Ez_[i][j][kGhostRight]=Ez_[i][j][kEnd_]; 

            Bx_[i][j][kGhostLeft]=Bx_[i][j][kBeg_]; 
            Bx_[i][j][kGhostRight]=Bx_[i][j][kEnd_]; 
            By_[i][j][kGhostLeft]=By_[i][j][kBeg_]; 
            By_[i][j][kGhostRight]=By_[i][j][kEnd_]; 
            Bz_[i][j][kGhostLeft]=Bz_[i][j][kBeg_]; 
            Bz_[i][j][kGhostRight]=Bz_[i][j][kEnd_]; 

            Jx_[i][j][kGhostLeft]=Jx_[i][j][kBeg_]; 
            Jx_[i][j][kGhostRight]=Jx_[i][j][kEnd_]; 
            Jy_[i][j][kGhostLeft]=Jy_[i][j][kBeg_]; 
            Jy_[i][j][kGhostRight]=Jy_[i][j][kEnd_]; 
            Jz_[i][j][kGhostLeft]=Jz_[i][j][kBeg_]; 
            Jz_[i][j][kGhostRight]=Jz_[i][j][kEnd_]; 
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
void Grid::getGhostVecAlt(const int side, double* ghostVec) {
    // this is the price of only 3 dimensions on field storage
    sliceMatToVec_(Ex_,side); 
    std::copy(sliceTmp_,sliceTmp_+nRealPtsYZPlane_,ghostVec+(0 * nRealPtsYZPlane_)); 

    sliceMatToVec_(Ey_,side); 
    std::copy(sliceTmp_,sliceTmp_+nRealPtsYZPlane_,ghostVec+(1 * nRealPtsYZPlane_)); 

    sliceMatToVec_(Ez_,side); 
    std::copy(sliceTmp_,sliceTmp_+nRealPtsYZPlane_,ghostVec+(2 * nRealPtsYZPlane_)); 

    sliceMatToVec_(Bx_,side); 
    std::copy(sliceTmp_,sliceTmp_+nRealPtsYZPlane_,ghostVec+(3 * nRealPtsYZPlane_)); 

    sliceMatToVec_(By_,side); 
    std::copy(sliceTmp_,sliceTmp_+nRealPtsYZPlane_,ghostVec+(4 * nRealPtsYZPlane_)); 

    sliceMatToVec_(Bz_,side); 
    std::copy(sliceTmp_,sliceTmp_+nRealPtsYZPlane_,ghostVec+(5 * nRealPtsYZPlane_)); 

    sliceMatToVec_(Jx_,side); 
    std::copy(sliceTmp_,sliceTmp_+nRealPtsYZPlane_,ghostVec+(6 * nRealPtsYZPlane_)); 

    sliceMatToVec_(Jy_,side); 
    std::copy(sliceTmp_,sliceTmp_+nRealPtsYZPlane_,ghostVec+(7 * nRealPtsYZPlane_)); 

    sliceMatToVec_(Jz_,side); 
    std::copy(sliceTmp_,sliceTmp_+nRealPtsYZPlane_,ghostVec+(8 * nRealPtsYZPlane_)); 
}; 

// alternate implementation of getGhostVec
// single loop, pulls from all fields at once
// worse performance due to more cache misses or negligible? 
void Grid::getGhostVec(const int side, double* ghostVec) {
    int i = sideToIndex_(side); 
    int j,k; // iterators
    int iter = -1; 
    for (j=jBeg_; j<jEnd_+1; ++j) { 
        for (k=kBeg_; k<kEnd_+1; ++k) { 
            ++iter; 
            ghostVec[0*nRealPtsYZPlane_ + iter]=Ex_[i][j][k];
            ghostVec[1*nRealPtsYZPlane_ + iter]=Ey_[i][j][k];
            ghostVec[2*nRealPtsYZPlane_ + iter]=Ez_[i][j][k];
            ghostVec[3*nRealPtsYZPlane_ + iter]=Bx_[i][j][k];
            ghostVec[4*nRealPtsYZPlane_ + iter]=By_[i][j][k];
            ghostVec[5*nRealPtsYZPlane_ + iter]=Bz_[i][j][k];
            ghostVec[6*nRealPtsYZPlane_ + iter]=Jx_[i][j][k];
            ghostVec[7*nRealPtsYZPlane_ + iter]=Jy_[i][j][k];
            ghostVec[8*nRealPtsYZPlane_ + iter]=Jz_[i][j][k];
        } 
    } 
};

// unbundles the data and puts it in the field 
// one loop for each field, aims to minimize cache misses
void Grid::setGhostVecAlt(const int side, const double* ghostVec) { 
    // this is the price of only 3 dimensions on field storage
    std::copy(ghostVec+(0 * nRealPtsYZPlane_),ghostVec+1*nRealPtsYZPlane_ -1,sliceTmp_); 
    unsliceMatToVec_(Ex_,side); 

    std::copy(ghostVec+(1 * nRealPtsYZPlane_),ghostVec+2*nRealPtsYZPlane_ -1,sliceTmp_); 
    unsliceMatToVec_(Ex_,side); 

    std::copy(ghostVec+(2 * nRealPtsYZPlane_),ghostVec+3*nRealPtsYZPlane_ -1,sliceTmp_); 
    unsliceMatToVec_(Ex_,side); 

    std::copy(ghostVec+(3 * nRealPtsYZPlane_),ghostVec+4*nRealPtsYZPlane_ -1,sliceTmp_); 
    unsliceMatToVec_(Ex_,side); 

    std::copy(ghostVec+(4 * nRealPtsYZPlane_),ghostVec+5*nRealPtsYZPlane_ -1,sliceTmp_); 
    unsliceMatToVec_(Ex_,side); 

    std::copy(ghostVec+(5 * nRealPtsYZPlane_),ghostVec+6*nRealPtsYZPlane_ -1,sliceTmp_); 
    unsliceMatToVec_(Ex_,side); 

    std::copy(ghostVec+(6 * nRealPtsYZPlane_),ghostVec+7*nRealPtsYZPlane_ -1,sliceTmp_); 
    unsliceMatToVec_(Ex_,side); 

    std::copy(ghostVec+(7 * nRealPtsYZPlane_),ghostVec+8*nRealPtsYZPlane_ -1,sliceTmp_); 
    unsliceMatToVec_(Ex_,side); 

    std::copy(ghostVec+(8 * nRealPtsYZPlane_),ghostVec+9*nRealPtsYZPlane_ -1,sliceTmp_); 
    unsliceMatToVec_(Ex_,side); 

};

// alternate implementation of setGhostVec
// single loop, sets all fields at once
// worse performance due to more cache misses or negligible? 
void Grid::setGhostVec(const int side, const double* ghostVec) {
    int i = sideToIndex_(side);
    int j,k; // iterators
    int iter = -1; 
    for (j=jBeg_; j<jEnd_+1; ++j) { 
        for (k=kBeg_; k<kEnd_+1; ++k) { 
            ++iter; 
            Ex_[i][j][k]=ghostVec[0*nRealPtsYZPlane_ + iter];
            Ey_[i][j][k]=ghostVec[1*nRealPtsYZPlane_ + iter];
            Ez_[i][j][k]=ghostVec[2*nRealPtsYZPlane_ + iter];
            Bx_[i][j][k]=ghostVec[3*nRealPtsYZPlane_ + iter];
            By_[i][j][k]=ghostVec[4*nRealPtsYZPlane_ + iter];
            Bz_[i][j][k]=ghostVec[5*nRealPtsYZPlane_ + iter];
            Jx_[i][j][k]=ghostVec[6*nRealPtsYZPlane_ + iter];
            Jy_[i][j][k]=ghostVec[7*nRealPtsYZPlane_ + iter];
            Jz_[i][j][k]=ghostVec[8*nRealPtsYZPlane_ + iter];
        } 
    } 
};

// function to convert -/+ 1 left/right side indicator to index in x direction 
int Grid::sideToIndex_(const int side) { 
    int i; 
    if (side < 0) { 
        i = iBeg_; 
    } 
    else { 
        i = iEnd_; 
    } 
    return i; 
}; 
