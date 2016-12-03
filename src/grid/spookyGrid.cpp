#include <stdlib.h>
#include <assert.h>
#include <math.h> 
#include <grid.hpp> 

// slices a physical plane in the x direction (excludes ghosts) 
void Grid::sliceMatToVec(const double*** mat, const int side) { 
    int j,k; // iterators
    if (side < 0) { 
        i = 1; 
    } 
    else { // side = +1
        i = nx_-(1+nGhosts_); 
    } 
    int iter=-1;
    for (j=1; j<ny_ - nGhosts_; ++j) { 
        for (k=1; k<nz_ - nGhosts_; ++k) { 
            // equivalently: iter = j*nz_ + k
            sliceTmp_[++iter] = mat[i][j][k]; 
        } 
    } 
}

// unslices a physical plane in the x direction (excludes ghosts) 
void Grid::unsliceMatToVec(const double*** mat, const int side) { 
    int j,k; // iterators
    if (side < 0) { 
        i = 1; 
    } 
    else { // side = +1
        i = nx_-(1+nGhosts_); 
    }
    int iter = -1; 
    for (j=1; j<ny_-nGhosts_; ++j) { 
        for (k=1; k<nz_-nGhosts_; ++k) { 
            // equivalently: iter = j*nz_ + k
            mat[i][j][k] = sliceTmp_[++iter]; 
        } 
    } 
} 

// updates ghost cells after evolving the field on physical points
// does this take a performance hit due to cache hit/miss?
void Grid::updateGhostCells() { 
    int i,j,k; // iterators 

    /* indices marking beginning and end of physical grids in each dimension */
    int iBeg = nGhosts_; 
    int iEnd = nx_ - (nGhosts_ + 1); 
    int jBeg = nGhosts_; 
    int jEnd = ny_ - (nGhosts_ + 1); 
    int kBeg = nGhosts_; 
    int kEnd = nz_ - (nGhosts_ + 1); 

    // update ghost cells in x direction 
    // iterates over yz plane
    for (j=1; j<ny_-nGhosts_; ++j) { 
        for (k=1; k<nz_-nGhosts_; ++k) { 
            Ex_[iBeg-1][j][k]=Ex_[iBeg][j][k]; 
            Ex_[iEnd+1][j][k]=Ex_[iEnd][j][k]; 
            Ey_[iBeg-1][j][k]=Ey_[iBeg][j][k]; 
            Ey_[iEnd+1][j][k]=Ey_[iEnd][j][k]; 
            Ez_[iBeg-1][j][k]=Ez_[iBeg][j][k]; 
            Ez_[iEnd+1][j][k]=Ez_[iEnd][j][k]; 

            Bx_[iBeg-1][j][k]=Bx_[iBeg][j][k]; 
            Bx_[iEnd+1][j][k]=Bx_[iEnd][j][k]; 
            By_[iBeg-1][j][k]=By_[iBeg][j][k]; 
            By_[iEnd+1][j][k]=By_[iEnd][j][k]; 
            Bz_[iBeg-1][j][k]=Bz_[iBeg][j][k]; 
            Bz_[iEnd+1][j][k]=Bz_[iEnd][j][k]; 

            Jx_[iBeg-1][j][k]=Jx_[iBeg][j][k]; 
            Jx_[iEnd+1][j][k]=Jx_[iEnd][j][k]; 
            Jy_[iBeg-1][j][k]=Jy_[iBeg][j][k]; 
            Jy_[iEnd+1][j][k]=Jy_[iEnd][j][k]; 
            Jz_[iBeg-1][j][k]=Jz_[iBeg][j][k]; 
            Jz_[iEnd+1][j][k]=Jz_[iEnd][j][k]; 
        } 
    }

    // updates ghost cells in y direction 
    // iterates over xz plane
    for (i=1; i<nx_-nGhosts_; ++i) { 
        for (k=1; k<nz_-nGhosts_; ++k) { 
            Ex_[i][jBeg-1][k]=Ex_[i][jBeg][k]; 
            Ex_[i][jEnd+1][k]=Ex_[i][jEnd][k]; 
            Ey_[i][jBeg-1][k]=Ey_[i][jBeg][k]; 
            Ey_[i][jEnd+1][k]=Ey_[i][jEnd][k]; 
            Ez_[i][jBeg-1][k]=Ez_[i][jBeg][k]; 
            Ez_[i][jEnd+1][k]=Ez_[i][jEnd][k]; 

            Bx_[i][jBeg-1][k]=Bx_[i][jBeg][k]; 
            Bx_[i][jEnd+1][k]=Bx_[i][jEnd][k]; 
            By_[i][jBeg-1][k]=By_[i][jBeg][k]; 
            By_[i][jEnd+1][k]=By_[i][jEnd][k]; 
            Bz_[i][jBeg-1][k]=Bz_[i][jBeg][k]; 
            Bz_[i][jEnd+1][k]=Bz_[i][jEnd][k]; 

            Jx_[i][jBeg-1][k]=Jx_[i][jBeg][k]; 
            Jx_[i][jEnd+1][k]=Jx_[i][jEnd][k]; 
            Jy_[i][jBeg-1][k]=Jy_[i][jBeg][k]; 
            Jy_[i][jEnd+1][k]=Jy_[i][jEnd][k]; 
            Jz_[i][jBeg-1][k]=Jz_[i][jBeg][k]; 
            Jz_[i][jEnd+1][k]=Jz_[i][jEnd][k]; 
        } 
    }

    // updates ghost cells in z direction 
    // iterates over xy plane 
    for (i=1; i<nx_-nGhosts_; ++i) { 
        for (j=1; j<ny_-nGhosts_; ++j) { 
            Ex_[i][j][kBeg-1]=Ex_[i][j][kBeg]; 
            Ex_[i][j][kEnd+1]=Ex_[i][j][kEnd]; 
            Ey_[i][j][kBeg-1]=Ey_[i][j][kBeg]; 
            Ey_[i][j][kEnd+1]=Ey_[i][j][kEnd]; 
            Ez_[i][j][kBeg-1]=Ez_[i][j][kBeg]; 
            Ez_[i][j][kEnd+1]=Ez_[i][j][kEnd]; 

            Bx_[i][j][kBeg-1]=Bx_[i][j][kBeg]; 
            Bx_[i][j][kEnd+1]=Bx_[i][j][kEnd]; 
            By_[i][j][kBeg-1]=By_[i][j][kBeg]; 
            By_[i][j][kEnd+1]=By_[i][j][kEnd]; 
            Bz_[i][j][kBeg-1]=Bz_[i][j][kBeg]; 
            Bz_[i][j][kEnd+1]=Bz_[i][j][kEnd]; 

            Jx_[i][j][kBeg-1]=Jx_[i][j][kBeg]; 
            Jx_[i][j][kEnd+1]=Jx_[i][j][kEnd]; 
            Jy_[i][j][kBeg-1]=Jy_[i][j][kBeg]; 
            Jy_[i][j][kEnd+1]=Jy_[i][j][kEnd]; 
            Jz_[i][j][kBeg-1]=Jz_[i][j][kBeg]; 
            Jz_[i][j][kEnd+1]=Jz_[i][j][kEnd]; 
        } 
    }

};

// returns size of ghost cell data to send 
int Grid::getGhostVecSize() { 
    // must bundle 9 scalar fields
    // each of size of the physical zy plane
    return 9*(ny_-2*nGhosts_)*(nz_-2*nGhosts_)
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
    sliceMatToVec(Ex_,side); 
    std::copy(sliceTmp_,sliceTmp_+ncp,ghostVec+0*ncp); 

    sliceMatToVec(Ey_,side); 
    std::copy(sliceTmp_,sliceTmp_+ncp,ghostVec+1*ncp); 

    sliceMatToVec(Ez_,side); 
    std::copy(sliceTmp_,sliceTmp_+ncp,ghostVec+2*ncp); 

    sliceMatToVec(Bx_,side); 
    std::copy(sliceTmp_,sliceTmp_+ncp,ghostVec+3*ncp); 

    sliceMatToVec(By_,side); 
    std::copy(sliceTmp_,sliceTmp_+ncp,ghostVec+4*ncp); 

    sliceMatToVec(Bz_,side); 
    std::copy(sliceTmp_,sliceTmp_+ncp,ghostVec+5*ncp); 

    sliceMatToVec(Jx_,side); 
    std::copy(sliceTmp_,sliceTmp_+ncp,ghostVec+6*ncp); 

    sliceMatToVec(Jy_,side); 
    std::copy(sliceTmp_,sliceTmp_+ncp,ghostVec+7*ncp); 

    sliceMatToVec(Jz_,side); 
    std::copy(sliceTmp_,sliceTmp_+ncp,ghostVec+8*ncp); 
}; 

// alternate implementation of getGhostVec
// single loop, pulls from all fields at once
// worse performance due to more cache misses or negligible? 
void Grid::getGhostVecAlt(const int side, double* ghostVec) {
    int realPts = (ny_ - 2*nGhosts_)*(nz_ - 2*nGhosts_); 
    if (side < 0) { 
        i = 1; 
    } 
    else { // side = +1
        i = nx_-(1+nGhosts_); 
    } 
    int j,k; // iterators
    int iter = -1; 
    for (j=1; j<ny_ - nGhosts_; ++j) { 
        for (k=1; k<nz_ - nGhosts_; ++k) { 
            ++iter; 
            ghostVec[0*realPts + iter]=Ex[i,j,k];
            ghostVec[1*realPts + iter]=Ey[i,j,k];
            ghostVec[2*realPts + iter]=Ez[i,j,k];
            ghostVec[3*realPts + iter]=Bx[i,j,k];
            ghostVec[4*realPts + iter]=By[i,j,k];
            ghostVec[5*realPts + iter]=Bz[i,j,k];
            ghostVec[6*realPts + iter]=Jx[i,j,k];
            ghostVec[7*realPts + iter]=Jy[i,j,k];
            ghostVec[8*realPts + iter]=Jz[i,j,k];
        } 
    } 
}

// unbundles the data and puts it in the field 
// one loop for each field, aims to minimize cache misses
int Grid::setGhostVec(const int side, const double* ghostVec) { 
    int ncp = nx_-2*nGhosts_; 
    // this is the price of only 3 dimensions on field storage
    std::copy(ghostVec+0*ncp,ghostVec+1*ncp -1,sliceTmp_); 
    unsliceMatToVec(Ex_,side); 

    std::copy(ghostVec+1*ncp,ghostVec+2*ncp -1,sliceTmp_); 
    unsliceMatToVec(Ex_,side); 

    std::copy(ghostVec+2*ncp,ghostVec+3*ncp -1,sliceTmp_); 
    unsliceMatToVec(Ex_,side); 

    std::copy(ghostVec+3*ncp,ghostVec+4*ncp -1,sliceTmp_); 
    unsliceMatToVec(Ex_,side); 

    std::copy(ghostVec+4*ncp,ghostVec+5*ncp -1,sliceTmp_); 
    unsliceMatToVec(Ex_,side); 

    std::copy(ghostVec+5*ncp,ghostVec+6*ncp -1,sliceTmp_); 
    unsliceMatToVec(Ex_,side); 

    std::copy(ghostVec+6*ncp,ghostVec+7*ncp -1,sliceTmp_); 
    unsliceMatToVec(Ex_,side); 

    std::copy(ghostVec+7*ncp,ghostVec+8*ncp -1,sliceTmp_); 
    unsliceMatToVec(Ex_,side); 

    std::copy(ghostVec+8*ncp,ghostVec+9*ncp -1,sliceTmp_); 
    unsliceMatToVec(Ex_,side); 

};

// alternate implementation of setGhostVec
// single loop, sets all fields at once
// worse performance due to more cache misses or negligible? 
void Grid::setGhostVecAlt(const int side, double* ghostVec) {
    int realPts = (ny_ - 2*nGhosts_)*(nz_ - 2*nGhosts_); 
    if (side < 0) { 
        i = 1; 
    } 
    else { // side = +1
        i = nx_-(1+nGhosts_); 
    } 
    int j,k; // iterators
    int iter = -1; 
    for (j=1; j<ny_ - nGhosts_; ++j) { 
        for (k=1; k<nz_ - nGhosts_; ++k) { 
            ++iter; 
            Ex[i,j,k]=ghostVec[0*realPts + iter];
            Ey[i,j,k]=ghostVec[1*realPts + iter];
            Ez[i,j,k]=ghostVec[2*realPts + iter];
            Bx[i,j,k]=ghostVec[3*realPts + iter];
            By[i,j,k]=ghostVec[4*realPts + iter];
            Bz[i,j,k]=ghostVec[5*realPts + iter];
            Jx[i,j,k]=ghostVec[6*realPts + iter];
            Jy[i,j,k]=ghostVec[7*realPts + iter];
            Jz[i,j,k]=ghostVec[8*realPts + iter];
        } 
    } 
}

