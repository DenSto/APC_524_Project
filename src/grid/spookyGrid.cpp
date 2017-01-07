#include <stdio.h>
#include <stdlib.h>
#include <algorithm> 
#include <assert.h> 
#include <math.h> 
#include "grid.hpp"

/// bundles the data in the ghost cells to send
/*! side = -/+ 1 for left/right x direction, -/+ 2 for y, -/+ 3 for z \n
 * ghostVec is the vector to store the data in, which must be of length ghostVecSize_ (can be determined with getGhostVecSize) \n
 * sendID = -2 to get Jrho fields, -1 to get EB fields, or sendID = an individual field ID (e.g. ExID_) to get just that field (used for Poisson updating for example) \n
 * Gets the data of the E,B,J fields along the specified boundary plane from the 1D array ghostVec to be sent with a single MPI call. If sendID = -1 (as used in each time step update), stores in order: Ex,Ey,Ez,Bx,By,Bz. If sendID = -2, stores in order: Jx,Jy,Jz,rho. \n
 * ghostVec can (and should) be unpacked with setGhostVec function 
 */ 
void Grid::getGhostVec(const int side, double* ghostVec, int sendID) {
    assert(-3 < sendID && sendID < nFieldsTotal_); 
    
    // create a temporary vector to store slices in 
    int n = maxPointsInPlane_;
    double* tmpVec = sliceTmp_; 

    // offset = 0 to get from the first/last physical cells 
    int offset=0;

    // determine number of fields being sent 
    int nfields; 
    switch (sendID) { 
        case -2: nfields=4; break; 
        case -1: nfields=6; break; 
        default: nfields=1; break; 
    }
    
    // "loop" over all fields to package 
    int begdex; 
    double*** field; 
    int fieldID,ifield;
    for (ifield=0; ifield<nfields; ++ifield) { 
        begdex=ifield*n; 
        switch (sendID) { 
            case -2: // send J/rho 
                switch (ifield) { 
                    case 0: fieldID = JxID_; break; 
                    case 1: fieldID = JyID_; break; 
                    case 2: fieldID = JzID_; break; 
                    case 3: fieldID = rhoID_; break; 
                }; 
                break; 
            case -1: // send E/B
                switch (ifield) { 
                    case 0: fieldID = ExID_; break; 
                    case 1: fieldID = EyID_; break; 
                    case 2: fieldID = EzID_; break; 
                    case 3: fieldID = BxID_; break; 
                    case 4: fieldID = ByID_; break; 
                    case 5: fieldID = BzID_; break; 
                }
                break; 
            default: fieldID = sendID; break; // send individual field 
        }; 
        field = fieldPtr_[fieldID]; 
        // slice the given field 
        sliceMatToVec_(fieldID,side,offset,tmpVec); 
        // store the slice in ghostVec 
        std::copy(tmpVec,tmpVec + n ,ghostVec + begdex); 
    } 
}; 

/// unbundles the data in the ghost cells to send
/*! side = -/+ 1 for left/right x direction, -/+ 2 for y, -/+ 3 for z \n
 * ghostVec is the vector to read the data from, which must be of length ghostVecSize_ (can be determined with getGhostVecSize) \n
 * sendID = -2 to set Jrho fields, -1 to set EB fields, or sendID = an individual field ID (e.g. ExID_) to set just that field (used for Poisson updating for example) \n
 * Sets the data of the E,B,J fields along the specified boundary plane from the 1D array ghostVec to be received with a single MPI call. If sendID = -1 (as used in each time step update), fields are read and set in order: Ex,Ey,Ez,Bx,By,Bz. If sendID = -2, fields are read and set in order: Jx,Jy,Jz,rho. \n
 * op is a flag determining how the field will be set. op = 0 replaces the current values in the field's ghost points with the values in ghostVec. op = 1 sums the current values and the values in ghostVec. \n
 * ghostVec can (and should) be generated with getGhostVec function 
 */ 
void Grid::setGhostVec(const int side, double* ghostVec, int sendID, int op) {
    assert(-3 < sendID && sendID < nFieldsTotal_); 
    
    // create a temporary vector to store slices in 
    int n = maxPointsInPlane_;
    double* tmpVec = sliceTmp_; 
    
    // offset = +1 to set into the RHS ghost vectors
    // offset = -1 to set into the LHS ghost vectors 
    int offset = 1; 
    if (side < 0) { 
        offset = -1; 
    };

    // determine number of fields being sent 
    int nfields; 
    switch (sendID) { 
        case -2: nfields=4; break; 
        case -1: nfields=6; break; 
        default: nfields=1; break; 
    }
    
    // "loop" over all fields to unpackage 
    int begdex,enddex; 
    double*** field; 
    int fieldID,ifield;
    for (ifield=0; ifield<nfields; ++ifield) { 
        begdex=ifield*n; 
        enddex=(ifield+1)*n; 
        // store the relevant portion for ghostVec into tmpVec
        std::copy(ghostVec + begdex, ghostVec + enddex, tmpVec);
        switch (sendID) { 
            case -2: // send J/rho 
                switch (ifield) { 
                    case 0: fieldID = JxID_; break; 
                    case 1: fieldID = JyID_; break; 
                    case 2: fieldID = JzID_; break; 
                    case 3: fieldID = rhoID_; break; 
                }; 
                break; 
            case -1: // send E/B
                switch (ifield) { 
                    case 0: fieldID = ExID_; break; 
                    case 1: fieldID = EyID_; break; 
                    case 2: fieldID = EzID_; break; 
                    case 3: fieldID = BxID_; break; 
                    case 4: fieldID = ByID_; break; 
                    case 5: fieldID = BzID_; break; 
                }
                break; 
            default: fieldID = sendID; break; // send individual field 
        }; 
        field = fieldPtr_[fieldID]; 
        // unslice the given field 
        unsliceMatToVec_(fieldID,side,offset,tmpVec,op); 
    } 
}; 


/// returns size of ghost cell data to send
/*! sendID is an integer specifying which fields are intended to be packaged into the ghost vector. \n 
 * -2: for J/rho package, -1 for E/B package, fieldID for any individual field (e.g. ExID_) \n
 * It is of length equal to the number of fields being sent times the maximum number of total points in any plane, so that it will be large enough to send the maximum amount of data in a single plane of any of the fields. 
 */
int Grid::getGhostVecSize(const int sendID) {
    assert(sendID > -3 && sendID < nFieldsTotal_);
    switch (sendID) {
        case -2: return 4*maxPointsInPlane_; break; // J,rho
        case -1: return 6*maxPointsInPlane_; break; // E,B
        default: return maxPointsInPlane_; break; // single field
    }
};

/// function to convert -/+ 1 left/right side indicator to index in x direction (description out of date) 
/*! For use with ghost cell methods. side=-1 indicates operations on the left side of the domain, side=+1 indicates operations on the right side of the domain. This method converts side into the correct index i to reference ghost cells on that side of the domain. For instance, called by getGhostVec and setGhostVec. Generalizes to any number of ghost cells so long as iBeg_ and iEnd_ are initialized correctly. 
 */
/// function to convert (-/+)(1,2,3) side indicator into (left/right)(x,y,z) index of boundary physical data point 
/*! Helper function for public ghost cell methods which accept side indicator as argument. \n 
 * Side < 0 will return index of first physical point, side > 0 will return index of last physical point \n 
 * abs(side) == 1 returns value in x direction, 2 in y, 3 in z \n 
 * This function is necessary because different field types have a different number of physical grid points in each direction. \n
 * fieldID is a private fieldID such as ExID_
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
/*! mat is 3D array whose real (non-ghost) data on one side will be replaced by data in the 1D array vec. vec must be of size maxPointsInPlane_. side is an integer -/+ 1 to indicate the location on the left/right side in the x direction, -/+ 2 in y, -/+ 3 in z. offset is an integer offset from the first/last physical index determined by side (e.g. side=-1 and offset=0 gives the yz plane of the 1st physical grid points in x direction, whereas offset=-1 would have returned the adjacent ghost cells and offset = 3 would have returned the 4th physical yz plane from the left). op=0 replaces the values in mat with those in vec, op=1 adds the values in vec to thos in mat. sliceMatToVec_ is the inverse function. 
 */ 
void Grid::unsliceMatToVec_(const int fieldID, const int side, const int offset, double* vec, const int op) { 
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
                switch (op) { 
                    case 0: mat[dex][j][k]=vec[++iter]; break; 
                    case 1: mat[dex][j][k]+=vec[++iter]; break; 
                }
            }
        } 
    } 
    else if (abs(side)==2) { 
        assert(dex < nyTot_); 
        for (i=iBeg_; i<iEnd; ++i) { 
            for (k=kBeg_; k<kEnd; ++k) { 
                switch (op) { 
                    case 0: mat[i][dex][k]=vec[++iter]; break; 
                    case 1: mat[i][dex][k]+=vec[++iter]; break; 
                } 
            }
        } 
    } 
    else if (abs(side)==3) { 
        assert(dex < nzTot_); 
        for (i=iBeg_; i<iEnd; ++i) { 
            for (j=jBeg_; j<jEnd; ++j) { 
                switch (op) { 
                    case 0: mat[i][j][dex]=vec[++iter]; break; 
                    case 1: mat[i][j][dex]+=vec[++iter]; break; 
                } 
            }
        } 
    } 
};

/// updates E,B ghost cells in y/z directions with periodic boundary conditions 
/*! Makes 4 calls each to get/setGhostVec for EB fields all at once
 */ 
void Grid::updatePeriodicGhostCells() { 
    // create a temporary vector to store ghostVecs 
    double* tmpGhost = ghostTmp_;  

    int sendID = -1;
    int op=0; 
    int side; 
    for (side=-3; side<4; ++side) { 
        // to set periodic boundary conditions in y/z directions, simply get/set ghostVec for side=+/-2, +/-3
        if (abs(side)>1) { 
            getGhostVec(side,tmpGhost,sendID); 
            setGhostVec(-side,tmpGhost,sendID,op); 
        };
    }; 
}; 
