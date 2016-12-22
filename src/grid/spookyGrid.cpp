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
    // create a temporary vector to store slices in 
    int n = maxPointsInPlane_; 
    double tmpVec[n]; 

    // offset = 0 to get from the first/last physical cells 
    int offset=0;

    // "loop" over all fields to package 
    int begdex; 
    double**** fieldPtr; 
    int fieldID; 
    int i; 
    for (i=0; i<nFieldsToSend_; ++i) { 
        begdex=i*n; 
        switch (i) { 
            case 0: fieldPtr = &Ex_; fieldID = ExID_; break; 
            case 1: fieldPtr = &Ey_; fieldID = EyID_; break; 
            case 2: fieldPtr = &Ez_; fieldID = EzID_; break; 
            case 3: fieldPtr = &Bx_; fieldID = BxID_; break; 
            case 4: fieldPtr = &By_; fieldID = ByID_; break; 
            case 5: fieldPtr = &Bz_; fieldID = BzID_; break; 
            case 6: fieldPtr = &Jx_; fieldID = JxID_; break; 
            case 7: fieldPtr = &Jy_; fieldID = JyID_; break; 
            case 8: fieldPtr = &Jz_; fieldID = JzID_; break; 
        }; 
        // slice the given field 
        sliceMatToVec_(*fieldPtr,side,fieldID,offset,tmpVec); 
        // store the slice in ghostVec 
        std::copy(tmpVec,tmpVec + n ,ghostVec + begdex); 
    };
}; 

/// unbundles the data in the ghost cells that have been received
/*! side = -/+ 1 for left/right x direction, -/+ 2 for y, -/+ 3 for z \n
 * ghostVec is the vector to read the data from, which must be of length ghostVecSize_ (can be determined with getGhostVecSize() ) \n
 * Sets the data of the E,B,J fields along the specified boundary plane from the 1D array received from a single MPI call. Sets in order: Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz. \n
 * ghostVec can be generated with getGhostVec function 
 */ 
void Grid::setGhostVec(const int side, double* ghostVec) {
    // create a temporary vector to store slices in 
    int n = maxPointsInPlane_;
    double tmpVec[n]; 
    
    // offset = +1 to set into the RHS ghost vectors
    // offset = -1 to set into the LHS ghost vectors 
    int offset = 1; 
    if (side < 0) { 
        offset = -1; 
    };

    // "loop" over all fields to unpackage 
    int begdex,enddex; 
    double**** fieldPtr; 
    int fieldID; 
    int i; 
    for (i=0; i<nFieldsToSend_; ++i) { 
        begdex=i*n; 
        enddex=(i+1)*n; 
        // store the relevant portion fo ghostVec into tmpVec
        std::copy(ghostVec + begdex, ghostVec + enddex, tmpVec);
        switch (i) { 
            case 0: fieldPtr = &Ex_; fieldID = ExID_; break; 
            case 1: fieldPtr = &Ey_; fieldID = EyID_; break; 
            case 2: fieldPtr = &Ez_; fieldID = EzID_; break; 
            case 3: fieldPtr = &Bx_; fieldID = BxID_; break; 
            case 4: fieldPtr = &By_; fieldID = ByID_; break; 
            case 5: fieldPtr = &Bz_; fieldID = BzID_; break; 
            case 6: fieldPtr = &Jx_; fieldID = JxID_; break; 
            case 7: fieldPtr = &Jy_; fieldID = JyID_; break; 
            case 8: fieldPtr = &Jz_; fieldID = JzID_; break; 
        }; 
        // unslice the given field 
        unsliceMatToVec_(*fieldPtr,side,fieldID,offset,tmpVec); 
    };
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
/*     
    sliceMatToVec_(Ex_,side,ExID_,0,tmpVec); 
    std::copy(tmpVec,tmpVec + n ,ghostVec + (++ifield * n)); 
   */  
    printf("suckers"); 
}; 
