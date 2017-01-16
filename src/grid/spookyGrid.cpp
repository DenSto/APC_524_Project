#include <stdio.h>
#include <stdlib.h>
#include <algorithm> 
#include <assert.h> 
#include <math.h> 
#include "grid.hpp"

/// bundles the data in the ghost cells to send
/*! side = -/+ i for left/right i-th direction
 *             i=1: x, i=2: y, i=3:, z \n
 * ghostVec is the vector to store the data in, which must be of length ghostVecSize_ \n
 * sendID >= 0: get an individual field ID (e.g. ExID_) \n
 *        = -1: (as used in each time step update), stores in order: Ex,Ey,Ez,Bx,By,Bz.\n 
 *        = -2: stores orthogonal J's and Rho in order: Jj,Jk,rho.\n
 * ghostVec can (and should) be unpacked with setGhostVec function 
 */ 
void Grid::getGhostVec(const int side, double* ghostVec, int sendID) {
    assert(-3 < sendID && sendID < nFieldsTotal_); 
    
    int offset; // offset from physical boundary
    if(sendID==-2){// get Jj,Jk,Rho on right ghost
        assert(side>0);
        offset = 1; // right ghost
    }else{// get physical values
        offset = 0; // physical
    }
    if(debug>1)fprintf(stderr,"rank=%d:get on side=%d with offset=%d,for field=%d\n",
                                  rank_MPI,side,offset,sendID);

    // create a temporary vector to store slices in 
    int n = maxPointsInPlane_;
    double* tmpVec = sliceTmp_; 

    const int xside = 1; 
    // const int yside = 2; // not used 
    const int zside = 3; 

    // determine number of fields being sent 
    int nfields; 
    switch (sendID) { 
        case -2: nfields=3; break; // two orthogonal J's and rho
        case -1: nfields=6; break; // two orthogonal B's should be sufficient
        default: nfields=1; break; // individual field
    }
    
    // "loop" over all fields to package 
    int begdex; 
    double*** field; 
    int fieldID;
    int ifield; 
    for (ifield=0; ifield<nfields; ++ifield) { 
        begdex=ifield*n; 
        switch (sendID) { 
            case -2: // send J/rho 
                switch (ifield) { 
                    // determine the 2 orthogonal J components to send 
                    case 0: 
                        if (abs(side) == xside) { 
                            fieldID = JyID_; 
                        } else { 
                            fieldID = JxID_; 
                        } 
                        break; 
                    case 1: 
                        if (abs(side) == zside) { 
                            fieldID = JyID_; 
                        } else { 
                            fieldID = JzID_; 
                        } 
                        break; 
                    case 2: fieldID = rhoID_; break; 
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

        // slice the given field with appropriate offset 
        sliceMatToVec_(fieldID,side,offset,tmpVec); 
        // store the slice in ghostVec 
        std::copy(tmpVec,tmpVec + n ,ghostVec + begdex); 
    } 
}; 

/// unbundles the data in the ghost cells to send
/*! side = -/+ i for left/right i-th direction \n
 * ghostVec is the vector to read the data from, which must be of length ghostVecSize_ \n
 * sendID = -2 to set Jj,Jk, rho fields
 *          -1 to set EB fields 
 *          an individual field ID (e.g. ExID_) to set just that field \n
 * ghostVec can (and should) be generated with getGhostVec function 
 */ 
void Grid::setGhostVec(const int side, double* ghostVec, int sendID) {
    assert(-3 < sendID && sendID < nFieldsTotal_); 
    
    int offset; // offset from physical boundary
    int option; // for unslice,0:replace, 1:sum
    if(sendID==-2){// sum Rho and Jj,Jk to left physical
        assert(side<0); //left
        offset = 0; // physical
        option = 1; // sum
    }else{// set ghost values by replace
        offset = abs(side)/side; // ghost
        option = 0; // replace
    }
    if(debug>1)fprintf(stderr,"rank=%d:set side=%d, offset=%d, field=%d, option=%d\n",
                                  rank_MPI,side,offset,sendID,option);

    // create a temporary vector to store slices in 
    int n = maxPointsInPlane_;
    double* tmpVec = sliceTmp_; 
    
    const int xside = 1; 
    // const int yside = 2; // not used 
    const int zside = 3; 
    
    // determine number of fields being sent 
    int nfields; 
    switch (sendID) { 
        case -2: nfields=3; break; // two orthogonal J's and rho
        case -1: nfields=6; break; // two orthogonal B's should be sufficient
        default: nfields=1; break; // individual field
    }
    
    // "loop" over all fields to unpackage 
    int begdex,enddex; 
    double*** field; 
    int fieldID; 
    int ifield;
    for (ifield=0; ifield<nfields; ++ifield) { 
        begdex=ifield*n; 
        enddex=(ifield+1)*n; 
        // store the relevant portion for ghostVec into tmpVec
        std::copy(ghostVec + begdex, ghostVec + enddex, tmpVec);
        switch (sendID) { 
            case -2: // set J/rho 
                switch (ifield) { 
                    // determine the 2 orthogonal J components to set 
                    case 0: 
                        if (abs(side) == xside) { 
                            fieldID = JyID_; 
                        } else { 
                            fieldID = JxID_; 
                        } 
                        break; 
                    case 1: 
                        if (abs(side) == zside) { 
                            fieldID = JyID_; 
                        } else { 
                            fieldID = JzID_; 
                        } 
                        break; 
                    case 2: fieldID = rhoID_; break; 
                }; 
                break; 
            case -1: // set E/B
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
            
        // unslice the given field with appropriate offset 
        unsliceMatToVec_(fieldID,side,offset,tmpVec,option);  
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
        case -2: return 3*maxPointsInPlane_; break; // rho, perp comps of J
        case -1: return 6*maxPointsInPlane_; break; // E,B
        default: return maxPointsInPlane_; break; // single field
    }
};

/// function to convert (-/+)(1,2,3) side indicator into (left/right)(x,y,z) index of boundary physical data point 
/*! Helper function for public ghost cell methods which accept side indicator as argument. \n 
 * Side < 0 will return index of first physical point, side > 0 will return index of last physical point \n
 * Notice that the last plane of a physical cell is a ghost plane! \n 
 * abs(side) == 1 returns value in x direction, 2 in y, 3 in z \n 
 * This function is necessary because different field types have a different number of grid points on physical cells in each direction. \n
 * However, notice that the last grid point on faces of physical cells is in fact a ghost point.
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

/// Gets offset of ghost points
/*! Called by getGhostVec. 
 * Returns the index of the value you want to get with getGhostVec \n 
 * side is the usual -/+ 1,2,3 for -/+ x,y,z boundaries \n 
 * fieldID is a fieldID (e.g. ExID_ from Grid)
 */ 
/*int Grid::getGhostOffset_(const int side, const int fieldID) { 
    assert(side != 0 && abs(side) < 4); 
    assert(-1 < fieldID && fieldID < nFieldsTotal_); 
    
    const int type = fieldType_[fieldID]; 
    
    const int dir = abs(side); 
    const int xdir = 1; 
    const int ydir = 2; 
    const int zdir = 3; 

    int offset = 0; 

    // get the magnitude of the offset 
    if (type == edgeXID_ && dir != xdir ) { 
        ++offset; 
    } else if (type == edgeYID_ && dir != ydir) { 
        ++offset; 
    } else if (type == edgeZID_ && dir != zdir) { 
        ++offset; 
    } else if (type == faceXID_ && dir == xdir) { 
        ++offset; 
    } else if (type == faceYID_ && dir == ydir) { 
        ++offset; 
    } else if (type == faceZID_ && dir == zdir) { 
        ++offset; 
    } else if (type == vertID_) { 
        ++offset; 
    } 
    
    // get the sign of the offset
    // opposite the sign of side 

    if (side > 1) offset *= -1; 
    return offset; 
};*/ 

/// slices a physical plane in the specified direction  
/*! mat is 3D array whose data on one side will be stored in vec as a 1D array. 
 * vec must be of size maxPointsInPlane_. 
 * side = 1: x right, -1: x left
 *      = 2: y right, -2: y left
 *      = 3: z right, -3: z left
 * offset is an integer offset from the first/last physical index 
 *      e.g. side=-1 and offset=0 gives the 1st physical yz plane in x direction
 *                       offset=-1 would have returned the adjacent ghost plane
 *                       offset = 3 would have returned the 4th physical yz plane 
 *           side=+1 and offset=0 gives the last physical yz plane in x direction.
 *
 * unsliceMatToVec_ is the inverse function. 
 */ 
void Grid::sliceMatToVec_(const int fieldID, const int side, const int offset, double* vec) { 
    // check for legal fieldID and side parameters 
    assert(fieldID > -1 && fieldID < nFieldsTotal_); 
    assert(side != 0 && abs(side) < 4); 

    // get the index to slice from
    int dex;
    if(side<0){//left side
        dex = nGhosts_; // index of first physical point, count from 0
    }else{
        dex = nGhosts_ + nxyzReal_[side-1] - 1; // index of last physical point
    }
    dex += offset;
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

/// unslices a plane in the specified direction  
/*! mat is 3D array on one sidei, will be replaced by data in the 1D array vec. 
 * vec must be of size maxPointsInPlane_. 
 * side = 1: x right, -1: x left
 *      = 2: y right, -2: y left
 *      = 3: z right, -3: z left
 * offset is an integer offset from the first/last physical index 
 *      e.g. side=-1 and offset=0 gives the 1st physical yz plane in x direction
 *                       offset=-1 would have returned the adjacent ghost plane
 *                       offset = 3 would have returned the 4th physical yz plane 
 *           side=+1 and offset=0 gives the last physical yz plane in x direction.
 * op = 0: replaces the values in mat with those in vec
 *      1: adds the values in vec to thos in mat. 
 *
 * sliceMatToVec_ is the inverse function. 
 */ 
void Grid::unsliceMatToVec_(const int fieldID, const int side, const int offset, double* vec, const int op) { 
    // check for legal fieldID and side parameters 
    assert(fieldID > -1 && fieldID < nFieldsTotal_); 
    assert(side != 0 && abs(side) < 4); 

    // get the index to unslice from
    int dex;
    if(side<0){//left side
        dex = nGhosts_; // index of first physical point, count from 0
    }else{
        dex = nGhosts_ + nxyzReal_[side-1] - 1; // index of last physical point
    }
    dex += offset;
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
/*! Makes 4 calls each to get/setGhostVec for EB fields all at once \n 
 * Deprecated method: correct performance not guaranteed 
 */ 
/*void Grid::updatePeriodicGhostCells() { 
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
};*/ 
