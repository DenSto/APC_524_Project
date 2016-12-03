#include <stdlib.h>
#include <assert.h>
#include <math.h> 
#include <grid.hpp> 

// slices a physical plane in the x direction (excludes ghosts) 
void sliceMatToVec(const double*** mat, const int side) { 
    int j,k; 
    if (side < 0) { 
        i = 1; 
    } 
    else { // side = +1
        i = nx_-(1+nGhosts_); 
    } 
    for (j=1; j<ny_ - nGhosts_; j++) { 
        for (k=1; k<nz_ - nGhosts_; k++) { 
            sliceTmp_[j*nz_ + k] = mat[i][j][k]; 
        } 
    } 
}

// unslices a physical plane in the x direction (excludes ghosts) 
void unsliceMatToVec(const double*** mat, const int side) { 
    int j,k; 
    if (side < 0) { 
        i = 1; 
    } 
    else { // side = +1
        i = nx_-(1+nGhosts_); 
    } 
    for (j=1; j<ny_-nGhosts_; j++) { 
        for (k=1; k<nz_-nGhosts_; k++) { 
            mat[i][j][k] = sliceTmp_[j*nz_ + k]; 
        } 
    } 
} 

// updates ghost cells after evolving the field on physical points
// not complete: needs to include all fields, and must resolve
// cache hit/miss design
void updateGhostCells() { 
    int i,j,k; 
    for (j=1; j<ny_-nGhosts_; j++) { 
        for (k=1; k<nz_-nGhosts_; k++) { 
            Ex_[0][j][k]=Ex_[1][j][k]; 
            Ex_[nx_-nGhosts_][j][k]=Ex_[1][j][k]; 
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
// returns error information (?) 
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

// unbundles the data and puts it in the field 
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
