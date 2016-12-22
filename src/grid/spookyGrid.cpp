#include <stdio.h>
#include <stdlib.h>
#include <algorithm> 
#include <assert.h> 
#include <math.h> 
#include "grid.hpp"

/// bundles the data in the ghost cells to send
/*! side = -/+ 1 for left/right x direction, -/+ 2 for y, -/+ 3 for z \n
 * ghostVec is the vector to store the data in, which must be of length ghostVecSize_ (can be determined with getGhostVecSize() ) \n
 * Stores the data of the E,B,J fields along the specified boundary plane into a 1D array to be sent with a single MPI call. Stores in order: Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz. \n
 * ghostVec can be unpacked with setGhostVec function 
 */ 
void Grid::getGhostVec(const int side, double* ghostVec) {
    int n = maxPointsInPlane_; 
    double tmpVec[n]; 
    int ifield=-1; 
    
    sliceMatToVec_(Ex_,side,ExID_,0,tmpVec); 
    std::copy(tmpVec,tmpVec + n ,ghostVec + (++ifield * n)); 

    sliceMatToVec_(Ey_,side,EyID_,0,tmpVec); 
    std::copy(tmpVec,tmpVec + n ,ghostVec + (++ifield * n)); 

    sliceMatToVec_(Ez_,side,EzID_,0,tmpVec); 
    std::copy(tmpVec,tmpVec + n ,ghostVec + (++ifield * n)); 

    sliceMatToVec_(Bx_,side,BxID_,0,tmpVec); 
    std::copy(tmpVec,tmpVec + n ,ghostVec + (++ifield * n)); 

    sliceMatToVec_(By_,side,ByID_,0,tmpVec); 
    std::copy(tmpVec,tmpVec + n ,ghostVec + (++ifield * n)); 

    sliceMatToVec_(Bz_,side,BzID_,0,tmpVec); 
    std::copy(tmpVec,tmpVec + n ,ghostVec + (++ifield * n)); 

    sliceMatToVec_(Jx_,side,JxID_,0,tmpVec); 
    std::copy(tmpVec,tmpVec + n ,ghostVec + (++ifield * n)); 

    sliceMatToVec_(Jy_,side,JyID_,0,tmpVec); 
    std::copy(tmpVec,tmpVec + n ,ghostVec + (++ifield * n)); 

    sliceMatToVec_(Jz_,side,JzID_,0,tmpVec); 
    std::copy(tmpVec,tmpVec + n ,ghostVec + (++ifield * n)); 
}; 

/// unbundles the data in the ghost cells that have been received
/*! side = -/+ 1 for left/right x direction, -/+ 2 for y, -/+ 3 for z \n
 * ghostVec is the vector to read the data from, which must be of length ghostVecSize_ (can be determined with getGhostVecSize() ) \n
 * Sets the data of the E,B,J fields along the specified boundary plane from the 1D array received from a single MPI call. Sets in order: Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz. \n
 * ghostVec can be generated with getGhostVec function 
 */ 
void Grid::setGhostVec(const int side, double* ghostVec) {
    // this is the price of only 3 dimensions on field storage
    int n = maxPointsInPlane_; 
    double tmpVec[n]; 
    int offset = 1; 
    if (side < 0) { 
        offset = -1; 
    }; 
    
    // commented out because "0" should be replaced with ifield
    // and "1" replaced with ++ifield, but compiler complained
    //int ifield=0; 
    
    std::copy(ghostVec + (0 * n), ghostVec + (1 * n), tmpVec); 
    unsliceMatToVec_(Ex_,side,ExID_,offset,tmpVec); 
    
    std::copy(ghostVec + (1 * n), ghostVec + (2 * n), tmpVec); 
    unsliceMatToVec_(Ey_,side,EyID_,offset,tmpVec); 
    
    std::copy(ghostVec + (2 * n), ghostVec + (3 * n), tmpVec); 
    unsliceMatToVec_(Ez_,side,EzID_,offset,tmpVec); 
    
    std::copy(ghostVec + (3 * n), ghostVec + (4 * n), tmpVec); 
    unsliceMatToVec_(Bx_,side,BxID_,offset,tmpVec); 
    
    std::copy(ghostVec + (4 * n), ghostVec + (5 * n), tmpVec); 
    unsliceMatToVec_(By_,side,ByID_,offset,tmpVec); 
    
    std::copy(ghostVec + (5 * n), ghostVec + (6 * n), tmpVec); 
    unsliceMatToVec_(Bz_,side,BzID_,offset,tmpVec); 
    
    std::copy(ghostVec + (6 * n), ghostVec + (7 * n), tmpVec); 
    unsliceMatToVec_(Jx_,side,JxID_,offset,tmpVec); 
    
    std::copy(ghostVec + (7 * n), ghostVec + (8 * n), tmpVec); 
    unsliceMatToVec_(Jy_,side,JyID_,offset,tmpVec); 
    
    std::copy(ghostVec + (8 * n), ghostVec + (9 * n), tmpVec); 
    unsliceMatToVec_(Jz_,side,JzID_,offset,tmpVec); 
}; 


/// returns size of ghost cell data to send
/*! this size is stored in the protected int ghostVecSize_
 */ 
int Grid::getGhostVecSize() { 
    // must bundle 9 scalar fields
    // each of size of the physical zy plane
    return ghostVecSize_; 
}; 

/// function to convert -/+ 1 left/right side indicator to index in x direction (description out of date) 
/*! For use with ghost cell methods. side=-1 indicates operations on the left side of the domain, side=+1 indicates operations on the right side of the domain. This method converts side into the correct index i to reference ghost cells on that side of the domain. For instance, called by getGhostVec and setGhostVec. Generalizes to any number of ghost cells so long as iBeg_ and iEnd_ are initialized correctly. 
 */ 
int Grid::sideToIndex_(const int side, const int fieldID) { 
    int dex; 
    if (side < 0) { 
        dex = 1; 
    }
    else { 
        dex = fieldSize_[fieldID][side-1]; 
    } 
    return dex; 
};

/// slices a physical plane in the specified direction (excludes ghosts) 
/*! mat is 3D array whose real (non-ghost) data on one side will be stored in vec as a 1D array. vec must be of size maxPointsInPlane_. side is an integer -/+ 1 to indicate the location on the left/right side in the x direction, -/+ 2 in y, -/+ 3 in z. offset is an integer offset from the first/last physical index determined by side (e.g. side=-1 and offset=0 gives the yz plane of the 1st physical grid points in x direction, whereas offset=-1 would have returned the adjacent ghost cells and offset = 3 would have returned the 4th physical yz plane from the left). unsliceMatToVec_ is the inverse function. 
 */ 
void Grid::sliceMatToVec_(double *** const mat, const int side, const int fieldID, const int offset, double* vec) { 
    assert(fieldID > -1 && fieldID < nIDs_); 
    assert(side != 0 && abs(side) < ndim_+1); 
    int dex = sideToIndex_(side,fieldID) + offset - (nGhosts_ + 1); 
    assert(dex > 0); 
    int i,j,k; // iterators
    int iter=-1; 

    if (abs(side)==1) {
        assert(dex < nxTot_); 
        int jEnd = fieldSize_[fieldID][1]; 
        int kEnd = fieldSize_[fieldID][2]; 
        for (j=jBeg_; j<jEnd; ++j) { 
            for (k=kBeg_; k<kEnd; ++k) { 
                vec[++iter] = mat[dex][j][k]; 
            }
        } 
    } 
    else if (abs(side)==2) { 
        assert(dex < nyTot_); 
        int iEnd = fieldSize_[fieldID][0]; 
        int kEnd = fieldSize_[fieldID][2]; 
        for (i=iBeg_; i<iEnd; ++i) { 
            for (k=kBeg_; k<kEnd; ++k) { 
                vec[++iter] = mat[i][dex][k]; 
            }
        } 
    } 
    else if (abs(side)==3) { 
        assert(dex < nzTot_); 
        int iEnd = fieldSize_[fieldID][0]; 
        int jEnd = fieldSize_[fieldID][1]; 
        for (i=iBeg_; i<iEnd; ++i) { 
            for (j=jBeg_; j<jEnd; ++j) { 
                vec[++iter] = mat[i][j][dex]; 
            }
        } 
    } 
};

/// unslices a physical plane in the specified direction (excludes ghosts) 
/*! mat is 3D array whose real (non-ghost) data on one side will be replaced by data in the 1D array vec. vec must be of size maxPointsInPlane_. side is an integer -/+ 1 to indicate the location on the left/right side in the x direction, -/+ 2 in y, -/+ 3 in z. offset is an integer offset from the first/last physical index determined by side (e.g. side=-1 and offset=0 gives the yz plane of the 1st physical grid points in x direction, whereas offset=-1 would have returned the adjacent ghost cells and offset = 3 would have returned the 4th physical yz plane from the left). sliceMatToVec_ is the inverse function. 
 */ 
void Grid::unsliceMatToVec_(double*** mat, const int side, const int fieldID, const int offset, double* vec) { 
    assert(fieldID > -1 && fieldID < nIDs_); 
    assert(side != 0 && abs(side) < ndim_+1); 
    int dex = sideToIndex_(side,fieldID) + offset - (nGhosts_ + 1); 
    assert(dex > 0); 
    int i,j,k; // iterators
    int iter=-1; 

    if (abs(side)==1) { 
        assert(dex < nxTot_); 
        int jEnd = fieldSize_[fieldID][1]; 
        int kEnd = fieldSize_[fieldID][2]; 
        for (j=jBeg_; j<jEnd; ++j) { 
            for (k=kBeg_; k<kEnd; ++k) { 
                mat[dex][j][k]=vec[++iter];
            }
        } 
    } 
    else if (abs(side)==2) { 
        assert(dex < nyTot_); 
        int iEnd = fieldSize_[fieldID][0]; 
        int kEnd = fieldSize_[fieldID][2]; 
        for (i=iBeg_; i<iEnd; ++i) { 
            for (k=kBeg_; k<kEnd; ++k) { 
                mat[i][dex][k]=vec[++iter]; 
            }
        } 
    } 
    else if (abs(side)==3) { 
        assert(dex < nzTot_); 
        int iEnd = fieldSize_[fieldID][0]; 
        int jEnd = fieldSize_[fieldID][1]; 
        for (i=iBeg_; i<iEnd; ++i) { 
            for (j=jBeg_; j<jEnd; ++j) { 
                mat[i][j][dex]=vec[++iter]; 
            }
        } 
    } 
};

// jlestz: fix this later to update all fields 
// making use of fieldSize_, sideToIndex_, etc. 
// (perhaps slice/unslicing)
void Grid::updatePeriodicGhostCells() { 
    printf("suckers"); 
}; 
