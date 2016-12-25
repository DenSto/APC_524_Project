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
// if sendID is -1, send JEB fields altogether 
// if sendID is any other value, send only that field
void Grid::getGhostVec(const int side, double* ghostVec, int sendID) {
    assert(sendID < nFieldsTotal_); 
    
    // create a temporary vector to store slices in 
    int n = maxPointsInPlane_;
    double* tmpVec = sliceTmp_; 

    // offset = 0 to get from the first/last physical cells 
    int offset=0;

    // "loop" over all fields to package 
    int begdex; 
    double*** field; 
    int fieldID; 

    if (sendID < 0) { 
    int i; 
        for (i=0; i<nFieldsToSend_; ++i) { 
            begdex=i*n; 
            switch (i) { 
                case 0: fieldID = ExID_; break; 
                case 1: fieldID = EyID_; break; 
                case 2: fieldID = EzID_; break; 
                case 3: fieldID = BxID_; break; 
                case 4: fieldID = ByID_; break; 
                case 5: fieldID = BzID_; break; 
                case 6: fieldID = JxID_; break; 
                case 7: fieldID = JyID_; break; 
                case 8: fieldID = JzID_; break; 
            }; 
            field = fieldPtr_[fieldID]; 
            // slice the given field 
            sliceMatToVec_(fieldID,side,offset,tmpVec); 
            // store the slice in ghostVec 
            std::copy(tmpVec,tmpVec + n ,ghostVec + begdex); 
        };
    } 
    else { 
        // slice the single field 
        sliceMatToVec_(fieldID,side,offset,tmpVec); 
        // store the slice in ghostVec 
        std::copy(tmpVec,tmpVec + n ,ghostVec); 
    } 
}; 

/// unbundles the data in the ghost cells that have been received
/*! side = -/+ 1 for left/right x direction, -/+ 2 for y, -/+ 3 for z \n
 * ghostVec is the vector to read the data from, which must be of length ghostVecSize_ (can be determined with getGhostVecSize() ) \n
 * Sets the data of the E,B,J fields along the specified boundary plane from the 1D array received from a single MPI call. Sets in order: Ex,Ey,Ez,Bx,By,Bz,Jx,Jy,Jz. \n
 * ghostVec can be generated with getGhostVec function 
 */ 
void Grid::setGhostVec(const int side, double* ghostVec, int sendID) {
    assert(sendID < nFieldsTotal_); 
    
    // create a temporary vector to store slices in 
    int n = maxPointsInPlane_;
    double* tmpVec = sliceTmp_; 
    
    // offset = +1 to set into the RHS ghost vectors
    // offset = -1 to set into the LHS ghost vectors 
    int offset = 1; 
    if (side < 0) { 
        offset = -1; 
    };

    // "loop" over all fields to unpackage 
    int begdex,enddex; 
    double*** field; 
    int fieldID; 
    
    if (sendID < 0) { 
        int i; 
        for (i=0; i<nFieldsToSend_; ++i) { 
            begdex=i*n; 
            enddex=(i+1)*n; 
            // store the relevant portion fo ghostVec into tmpVec
            std::copy(ghostVec + begdex, ghostVec + enddex, tmpVec);
            switch (i) { 
                case 0: fieldID = ExID_; break; 
                case 1: fieldID = EyID_; break; 
                case 2: fieldID = EzID_; break; 
                case 3: fieldID = BxID_; break; 
                case 4: fieldID = ByID_; break; 
                case 5: fieldID = BzID_; break; 
                case 6: fieldID = JxID_; break; 
                case 7: fieldID = JyID_; break; 
                case 8: fieldID = JzID_; break; 
            }; 
            field = fieldPtr_[fieldID]; 
            // unslice the given field 
            unsliceMatToVec_(fieldID,side,offset,tmpVec); 
        };
    } 
    else { 
        // store the slice in ghostVec 
        std::copy(ghostVec,ghostVec + n ,tmpVec); 
        // unslice the single field 
        unsliceMatToVec_(fieldID,side,offset,tmpVec); 
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

/// function to convert -/+ 1 left/right side indicator to index in x direction (description out of date) 
/*! For use with ghost cell methods. side=-1 indicates operations on the left side of the domain, side=+1 indicates operations on the right side of the domain. This method converts side into the correct index i to reference ghost cells on that side of the domain. For instance, called by getGhostVec and setGhostVec. Generalizes to any number of ghost cells so long as iBeg_ and iEnd_ are initialized correctly. 
 */ 
int Grid::sideToIndex_(const int side, const int fieldID) { 
    int dex; 
    if (side < 0) { 
        dex = 1; 
    }
    else { 
        int type = fieldType_[fieldID]; 
        // since side > 0 in this branch, side-1 converts (x,y,z = 1,2,3) --> (0,1,2)
        int dir = side-1;
        // fieldSize_ is the total number of points
        // -1 to convert to 0 indexing 
        // -nGhosts to subtract ghost cells 
        dex = fieldSize_[type][dir]-(nGhosts_+1); 
    } 
    return dex; 
};

/// slices a physical plane in the specified direction (excludes ghosts) 
/*! mat is 3D array whose real (non-ghost) data on one side will be stored in vec as a 1D array. vec must be of size maxPointsInPlane_. side is an integer -/+ 1 to indicate the location on the left/right side in the x direction, -/+ 2 in y, -/+ 3 in z. offset is an integer offset from the first/last physical index determined by side (e.g. side=-1 and offset=0 gives the yz plane of the 1st physical grid points in x direction, whereas offset=-1 would have returned the adjacent ghost cells and offset = 3 would have returned the 4th physical yz plane from the left). unsliceMatToVec_ is the inverse function. 
 */ 
void Grid::sliceMatToVec_(const int fieldID, const int side, const int offset, double* vec) { 
    // check for legal fieldID and side parameters 
    assert(fieldID > -1 && fieldID < nFieldsTotal_); 
    assert(side != 0 && abs(side) < 4); 

    // get the index to slice from
    int dex = sideToIndex_(side,fieldID) + offset; 
    assert(dex > -1); 
    
    // use fieldID to get the pointer to the field 
    double*** mat = fieldPtr_[fieldID]; 
   
    // directions 
    int xside=1; 
    int yside=2; 
    int zside=3; 
    
    // limits 
    int iEnd = sideToIndex_(xside,fieldID)+1;  
    int jEnd = sideToIndex_(yside,fieldID)+1; 
    int kEnd = sideToIndex_(zside,fieldID)+1; 
    
    // iterators 
    int i,j,k; 
    int iter=-1;

    if (abs(side)==1) {
        assert(dex < nxTot_);
        for (j=jBeg_; j<jEnd; ++j) { 
            for (k=kBeg_; k<kEnd; ++k) { 
                vec[++iter] = mat[dex][j][k]; 
            }
        } 
    } 
    else if (abs(side)==2) { 
        assert(dex < nyTot_); 
        for (i=iBeg_; i<iEnd; ++i) { 
            for (k=kBeg_; k<kEnd; ++k) { 
                vec[++iter] = mat[i][dex][k]; 
            }
        } 
    } 
    else if (abs(side)==3) { 
        assert(dex < nzTot_); 
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
void Grid::unsliceMatToVec_(const int fieldID, const int side, const int offset, double* vec) { 
    // check for legal fieldID and side parameters 
    assert(fieldID > -1 && fieldID < nFieldsTotal_); 
    assert(side != 0 && abs(side) < 4); 

    // get the index to unslice from
    int dex = sideToIndex_(side,fieldID) + offset; 
    assert(dex > -1); 

    // use fieldID to get the pointer to the field 
    double*** mat = fieldPtr_[fieldID]; 
    
    // directions 
    int xside=1; 
    int yside=2; 
    int zside=3; 
    
    // limits 
    int iEnd = sideToIndex_(xside,fieldID)+1;  
    int jEnd = sideToIndex_(yside,fieldID)+1; 
    int kEnd = sideToIndex_(zside,fieldID)+1; 
    
    // iterators 
    int i,j,k; 
    int iter=-1; 

    if (abs(side)==1) { 
        assert(dex < nxTot_); 
        for (j=jBeg_; j<jEnd; ++j) { 
            for (k=kBeg_; k<kEnd; ++k) { 
                mat[dex][j][k]=vec[++iter];
            }
        } 
    } 
    else if (abs(side)==2) { 
        assert(dex < nyTot_); 
        for (i=iBeg_; i<iEnd; ++i) { 
            for (k=kBeg_; k<kEnd; ++k) { 
                mat[i][dex][k]=vec[++iter]; 
            }
        } 
    } 
    else if (abs(side)==3) { 
        assert(dex < nzTot_); 
        for (i=iBeg_; i<iEnd; ++i) { 
            for (j=jBeg_; j<jEnd; ++j) { 
                mat[i][j][dex]=vec[++iter]; 
            }
        } 
    } 
};

/// updates J,E,B ghost cells in y/z directions with periodic boundary conditions 
void Grid::updatePeriodicGhostCells() { 
    // create a temporary vector to store ghostVecs 
    double* tmpGhost = ghostTmp_;  

    int sendID = -1; 
    int side; 
    for (side=-3; side<4; ++side) { 
        // to set periodic boundary conditions in y/z directions, simply get/set ghostVec for side=+/-2, +/-3
        if (abs(side)>1) { 
            getGhostVec(side,tmpGhost,sendID); 
            setGhostVec(-side,tmpGhost,sendID); 
        };
    }; 
}; 
