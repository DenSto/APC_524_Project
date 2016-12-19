#include <stdio.h>
#include <stdlib.h>
#include <algorithm> // used for std::copy
#include <assert.h> // plans to use?
#include <math.h> // unneeded?
#include "grid.hpp"

/// slices a physical plane in the x direction (excludes ghosts) 
/*! mat is 3D array whose real (non-ghost) data on one side will be stored in sliceTmp_ as a 1D array. side is an integer +1 to indicate storage of the right hand side values and -1 to indicate storage of the left hand side. Complementary function to unsliceMatToVec_. 
 */ 
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

/// unslices a physical plane in the x direction (excludes ghosts) 
/*! mat is 3D array whose real (non-ghost) data on one side will be set from the temporary 1D array sliceTmp_. side is an integer +1 to indicate setting of the right hand side values and -1 to indicate setting of the left hand side. Complementary function to sliceMatToVec_. 
 */ 
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

/// updates ghost cells after evolving the field on physical points
/*! For each of Ei_,Bi_,Ji_ (for i=x,y,z), updates the ghost cells in the y and z directions with periodic boundary conditions. \n
 * Updates of ghost cells in the x direction requires MPI calls due to 1D domain decomposition, and is handled in domain class, not here. \n
 * Currently this method requires nGhosts_=1 and will not perform correctly if nGhosts_ != 1 (it may not crash but will not update the ghost cells as desired). 
 */ 
void Grid::updatePeriodicGhostCells() { 
    int i,j,k; // iterators 
    int jGhostLeft=jBeg_-1; 
    int jGhostRight=jEnd_+1; 
    int kGhostLeft=kBeg_-1; 
    int kGhostRight=kEnd_+1; 
    
    // updates ghost cells in y direction 
    // iterates over xz plane
    for (i=iBeg_; i<iEnd_+1; ++i) { 
        for (k=kBeg_; k<kEnd_+1; ++k) { 
            Ex_[i][jGhostLeft][k]=Ex_[i][jEnd_][k]; 
            Ex_[i][jGhostRight][k]=Ex_[i][jBeg_][k]; 
            Ey_[i][jGhostLeft][k]=Ey_[i][jEnd_][k]; 
            Ey_[i][jGhostRight][k]=Ey_[i][jBeg_][k]; 
            Ez_[i][jGhostLeft][k]=Ez_[i][jEnd_][k]; 
            Ez_[i][jGhostRight][k]=Ez_[i][jBeg_][k]; 

            Bx_[i][jGhostLeft][k]=Bx_[i][jEnd_][k]; 
            Bx_[i][jGhostRight][k]=Bx_[i][jBeg_][k]; 
            By_[i][jGhostLeft][k]=By_[i][jEnd_][k]; 
            By_[i][jGhostRight][k]=By_[i][jBeg_][k]; 
            Bz_[i][jGhostLeft][k]=Bz_[i][jEnd_][k]; 
            Bz_[i][jGhostRight][k]=Bz_[i][jBeg_][k]; 

            Jx_[i][jGhostLeft][k]=Jx_[i][jEnd_][k]; 
            Jx_[i][jGhostRight][k]=Jx_[i][jBeg_][k]; 
            Jy_[i][jGhostLeft][k]=Jy_[i][jEnd_][k]; 
            Jy_[i][jGhostRight][k]=Jy_[i][jBeg_][k]; 
            Jz_[i][jGhostLeft][k]=Jz_[i][jEnd_][k]; 
            Jz_[i][jGhostRight][k]=Jz_[i][jBeg_][k]; 
        } 
    }

    // updates ghost cells in z direction 
    // iterates over xy plane 
    for (i=iBeg_; i<iEnd_+1; ++i) { 
        for (j=jBeg_; j<jEnd_+1; ++j) { 
            Ex_[i][j][kGhostLeft]=Ex_[i][j][kEnd_]; 
            Ex_[i][j][kGhostRight]=Ex_[i][j][kBeg_]; 
            Ey_[i][j][kGhostLeft]=Ey_[i][j][kEnd_]; 
            Ey_[i][j][kGhostRight]=Ey_[i][j][kBeg_]; 
            Ez_[i][j][kGhostLeft]=Ez_[i][j][kEnd_]; 
            Ez_[i][j][kGhostRight]=Ez_[i][j][kBeg_]; 

            Bx_[i][j][kGhostLeft]=Bx_[i][j][kEnd_]; 
            Bx_[i][j][kGhostRight]=Bx_[i][j][kBeg_]; 
            By_[i][j][kGhostLeft]=By_[i][j][kEnd_]; 
            By_[i][j][kGhostRight]=By_[i][j][kBeg_]; 
            Bz_[i][j][kGhostLeft]=Bz_[i][j][kEnd_]; 
            Bz_[i][j][kGhostRight]=Bz_[i][j][kBeg_]; 

            Jx_[i][j][kGhostLeft]=Jx_[i][j][kEnd_]; 
            Jx_[i][j][kGhostRight]=Jx_[i][j][kBeg_]; 
            Jy_[i][j][kGhostLeft]=Jy_[i][j][kEnd_]; 
            Jy_[i][j][kGhostRight]=Jy_[i][j][kBeg_]; 
            Jz_[i][j][kGhostLeft]=Jz_[i][j][kEnd_]; 
            Jz_[i][j][kGhostRight]=Jz_[i][j][kBeg_]; 
        } 
    }

};

/// returns size of ghost cell data to send
/*! this size is stored in the protected int ghostVecSize_
 */ 
int Grid::getGhostVecSize() { 
    // must bundle 9 scalar fields
    // each of size of the physical zy plane
    return ghostVecSize_; 
}; 

/// bundles the data in the ghost cells to send
/*! stores the data of the E,B,J fields at all of the ghost points along the domain interfaces (yz plane) into a 1D array of doubles to be sent with a single MPI call. ghostVec is an array of length ghostVecSize_ to store the data in. Side is -1 for left side of domain, +1 for right side. Sends (in order): Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz. This is an alternative implementation to the one in getGhostVecAlt which is less elegant but might decrease cache misses? Requires profiling
 */ 
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

/// bundles the data in the ghost cells to send
/*! stores the data of the E,B,J fields at all of the ghost points along the domain interfaces (yz plane) into a 1D array of doubles to be sent with a single MPI call. ghostVec is an array of length ghostVecSize_ to store the data in. Side is -1 for left side of domain, +1 for right side. Sends (in order): Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz. This is a more elegant implementation than the one in getGhostVec, but may increase cache misses? Requires profiling. 
 */ 
void Grid::getGhostVec(const int side, double* ghostVec) {
    int i = sideToIndex_(side); 
    int j,k; // iterators
    int iter = -1; 
    for (j=jBeg_; j<jEnd_+1; ++j) { 
        for (k=kBeg_; k<kEnd_+1; ++k) { 
            ++iter;
            //fprintf(stderr,"iter=%d\n",iter); 
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

///  unbundles the data sent from the ghost cells and puts it in the field 
/*! to be used in conjuction with getGhostVec or getGhostVecAlt. ghostVec is a 1D array storing each of the E,B,J field values at each of the ghost points along the domain interfaces (yz plane) of a single side. Side specifies which side this data should be set to. -1 corresponds to the left side of the domain (receiving data from the right ghost cells of the previous domain) and +1 to the right side (receiving data from the left ghost cells of the next domain). This is an alternate implementation of setGhostVec. setGhostVecAlt is less elegant but may reduce cache misses (requires profiling). 
 */ 
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

///  unbundles the data sent from the ghost cells and puts it in the field 
/*! to be used in conjuction with getGhostVec or getGhostVecAlt. ghostVec is a 1D array storing each of the E,B,J field values at each of the ghost points along the domain interfaces (yz plane) of a single side. Side specifies which side this data should be set to. -1 corresponds to the left side of the domain (receiving data from the right ghost cells of the previous domain) and +1 to the right side (receiving data from the left ghost cells of the next domain). This is an alternate implementation of setGhostVecAlt. setGhostVec is more elegant but may increase cache misses (requires profiling). 
 */ 
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

/// function to convert -/+ 1 left/right side indicator to index in x direction 
/*! For use with ghost cell methods. side=-1 indicates operations on the left side of the domain, side=+1 indicates operations on the right side of the domain. This method converts side into the correct index i to reference ghost cells on that side of the domain. For instance, called by getGhostVec and setGhostVec. Generalizes to any number of ghost cells so long as iBeg_ and iEnd_ are initialized correctly. 
 */ 
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
