#define MAIN_CPP
#include "gtest/gtest.h"
#include "grid.hpp"
#include "math.h"
#include "RNG.hpp"

// Tests of input/output to public methods
class GridTest : public ::testing::Test {
   protected:
      virtual void SetUp() {
         int nxyz [3] = {4,6,8};  
         int nGhosts = 1; 
         double xyz0 [3] = {0,0,0};
         double Lxyz [3] = {1,1,1};

         grid = new Grid(nxyz, nGhosts, xyz0, Lxyz);

      }

      virtual void TearDown() {
         delete grid;
      }

      Grid *grid;

};

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}


// tests requiring internal examination
class GridPrivateTest : public ::testing::Test {
   protected:
      virtual void SetUp() {
         int nxyz [3] = {4,6,8};
         int nGhosts = 1; 
         double xyz0 [3] = {0,0,0};
         double Lxyz [3] = {1,1,1};

         grid = new Grid(nxyz, nGhosts, xyz0, Lxyz);

      }

      virtual void TearDown() {
         delete grid;
      }

      double sumField(double*** field) { 
          int i,j,k;
          double fieldSum = 0; 
          for (i=0; i<grid->nxTot_; ++i) { 
              for (j=0; j<grid->nyTot_; ++j) { 
                  for (k=0; k<grid->nzTot_; ++k) { 
                      fieldSum += field[i][j][k]; 
                  }
              }
          }
          return fieldSum; 
      }; 

      Grid *grid;

};

// test sideToIndex
TEST_F(GridPrivateTest, sideToIndexTest) { 
    int ID = grid->ExID_; 

    int xside=1; 
    int yside=2; 
    int zside=3; 

    // hard-coded values for Ex
    int xLeftReal = 1; 
    int yLeftReal = 1; 
    int zLeftReal = 1; 
    int xRightReal = grid->nx_-2; 
    int yRightReal = grid->ny_-1; 
    int zRightReal = grid->nz_-1; 

    EXPECT_EQ(xLeftReal,grid->sideToIndex_(-xside,ID)); 
    EXPECT_EQ(xRightReal,grid->sideToIndex_(xside,ID)); 
    EXPECT_EQ(yLeftReal,grid->sideToIndex_(-yside,ID)); 
    EXPECT_EQ(yRightReal,grid->sideToIndex_(yside,ID)); 
    EXPECT_EQ(zLeftReal,grid->sideToIndex_(-zside,ID)); 
    EXPECT_EQ(zRightReal,grid->sideToIndex_(zside,ID)); 
}; 

// tests updateGhostCells for one of the fields 
// passing this test also implies that slice/unslice methods work
// passing this test also implies that get/set ghost methods work 
TEST_F(GridPrivateTest, periodicUpdateTest) {
    // choose to test Jz (all of JEB are sent in updatePeriodic)
    int ID = grid->ExID_; 
    double*** field = grid->fieldPtr_[ID]; 

    // store incremented value at each point 
    int i,j,k; 
    double iter=0; 
    for (i = 0; i<grid->nxTot_; ++i) { 
        for (j=0; j<grid->nyTot_; ++j) { 
            for (k=0; k<grid->nzTot_; ++k) { 
                ++iter; 
                field[i][j][k] += pow(iter,2); 
            } 
        } 
    } 

    // debug: verify the field is being filled up as expected 
    // this only applies if field is being loaded with += iteri
    // since this is an analytic sum of the total number of elements
    /* 
    int n = grid->nxTot_*grid->nyTot_*grid->nzTot_; 
    int anSum = n*(n+1)/2; 
    double sumTot = sumField(field); 
    EXPECT_EQ(sumTot,anSum); 
    */ 

    // get the initial sums of the relevant rows and columns
    double sumYLeftGhost = 0; 
    double sumYLeftReal = 0; 
    double sumYRightGhost = 0; 
    double sumYRightReal = 0;

    double sumZLeftGhost = 0; 
    double sumZLeftReal = 0; 
    double sumZRightGhost = 0; 
    double sumZRightReal = 0; 

    // indices for directions 
    int xside=1; 
    int yside=2; 
    int zside=3; 

    int dexXLeftReal = grid->sideToIndex_(-xside,ID); 
    int dexXRightReal = grid->sideToIndex_(xside,ID); 
    
    int dexYLeftReal = grid->sideToIndex_(-yside,ID); 
    int dexYRightReal = grid->sideToIndex_(yside,ID); 
    int dexYLeftGhost = dexYLeftReal - 1; 
    int dexYRightGhost = dexYRightReal + 1; 

    int dexZLeftReal = grid->sideToIndex_(-zside,ID); 
    int dexZRightReal = grid->sideToIndex_(zside,ID); 
    int dexZLeftGhost = dexZLeftReal - 1; 
    int dexZRightGhost = dexZRightReal + 1; 

    double val = 0;
    int jreal, kreal; 
    for (i = dexXLeftReal; i < dexXRightReal+1; ++i) { 
        for (j = dexYLeftGhost; j < dexYRightGhost+1; ++j) { 
            if (j==dexYLeftGhost || j==dexYRightGhost) { 
                jreal=0; 
            } 
            else { 
                jreal=1;
            } 
            for (k = dexZLeftGhost; k < dexZRightGhost+1; ++k) { 
                if (k==dexZLeftGhost || k==dexZRightGhost) { 
                    kreal=0; 
                } 
                else { 
                    kreal=1;
                } 

                val = field[i][j][k]; 
                
                // sums for constant y
                if (kreal==1) { 
                    if (j==dexYLeftGhost) { 
                        sumYLeftGhost += val; 
                    } 
                    else if (j==dexYLeftReal) { 
                        sumYLeftReal += val; 
                    } 
                    else if (j==dexYRightReal) { 
                        sumYRightReal += val; 
                    } 
                    else if (j==dexYRightGhost) { 
                        sumYRightGhost += val; 
                    }
                }

                // sums for constant z 
                if (jreal==1) { 
                    if (k==dexZLeftGhost) { 
                        sumZLeftGhost += val; 
                    } 
                    else if (k==dexZLeftReal) { 
                        sumZLeftReal += val; 
                    } 
                    else if (k==dexZRightReal) { 
                        sumZRightReal += val; 
                    } 
                    else if (k==dexZRightGhost) { 
                        sumZRightGhost += val; 
                    }
                }
            } 
        } 
    }

    // impose periodic BC in y and z 
    grid->updatePeriodicGhostCells(); 

    // sums after periodic boundary conditions 
    double newSumYLeftGhost = 0; 
    double newSumYLeftReal = 0; 
    double newSumYRightGhost = 0; 
    double newSumYRightReal = 0;

    double newSumZLeftGhost = 0; 
    double newSumZLeftReal = 0; 
    double newSumZRightGhost = 0; 
    double newSumZRightReal = 0; 

    for (i = dexXLeftReal; i < dexXRightReal+1; ++i) { 
        for (j = dexYLeftGhost; j < dexYRightGhost+1; ++j) { 
            if (j==dexYLeftGhost || j==dexYRightGhost) { 
                jreal=0; 
            } 
            else { 
                jreal=1;
            } 
            for (k = dexZLeftGhost; k < dexZRightGhost+1; ++k) { 
                if (k==dexZLeftGhost || k==dexZRightGhost) { 
                    kreal=0; 
                } 
                else { 
                    kreal=1;
                } 
                
                val = field[i][j][k]; 
                
                // newSums for constant y 
                if (kreal==1) { 
                    if (j==dexYLeftGhost) { 
                        newSumYLeftGhost += val; 
                    } 
                    else if (j==dexYLeftReal) { 
                        newSumYLeftReal += val; 
                    } 
                    else if (j==dexYRightReal) { 
                        newSumYRightReal += val; 
                    } 
                    else if (j==dexYRightGhost) { 
                        newSumYRightGhost += val; 
                    }
                }

                // newSums for constant z 
                if (jreal==1) { 
                    if (k==dexZLeftGhost) { 
                        newSumZLeftGhost += val; 
                    } 
                    else if (k==dexZLeftReal) { 
                        newSumZLeftReal += val; 
                    } 
                    else if (k==dexZRightReal) { 
                        newSumZRightReal += val; 
                    } 
                    else if (k==dexZRightGhost) { 
                        newSumZRightGhost += val; 
                    }
                }
            } 
        } 
    }

    // check that BCs were updated appropriately 
    EXPECT_EQ(newSumYLeftGhost,sumYRightReal); 
    EXPECT_EQ(newSumYRightGhost,sumYLeftReal); 
    EXPECT_EQ(newSumZLeftGhost,sumZRightReal); 
    EXPECT_EQ(newSumZRightGhost,sumZLeftReal); 

    // check that real cells were unaffected 
    EXPECT_EQ(newSumYLeftReal,sumYLeftReal); 
    EXPECT_EQ(newSumYRightReal,sumYRightReal); 
    EXPECT_EQ(newSumZLeftReal,sumZLeftReal); 
    EXPECT_EQ(newSumZRightReal,sumZRightReal); 
}

// test correct ghostVecSize
// compile error: says nxyz is out of scope even though it is in SetUp()?
TEST_F(GridPrivateTest, ghostVecSizeTest) { 
    int nx = grid->nxTot_; 
    int ny = grid->nyTot_; 
    int nz = grid->nzTot_; 

    int maxPlane = std::max(nx*ny,ny*nz); 
    maxPlane = std::max(maxPlane,nx*nz); 

    EXPECT_EQ(maxPlane,grid->maxPointsInPlane_);

    int nFields = 6; 
    int nSources= 4; 

    int sizeF = nFields*maxPlane; 
    int sizeS = nSources*maxPlane; 

    EXPECT_EQ(sizeF, grid->getGhostVecSize(-1)); 
    EXPECT_EQ(sizeS, grid->getGhostVecSize(-2)); 
} 


