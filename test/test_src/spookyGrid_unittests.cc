#include "gtest/gtest.h"
#include "grid.hpp"
#include "math.h"

// Tests of input/output to public methods
class GridTest : public ::testing::Test {
   protected:
      virtual void SetUp() {
         int nxyz [3] = {3,5,7};
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
class gridPrivateTest : public ::testing::Test {
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

// test updatePeriodicGhostCells in Ex
TEST_F(gridPrivateTest, periodicUpdateTest) {
    // choose to test Ex
    int ID = grid->ExID_; 
    int type = grid->fieldType_[ID]; 
    double*** field = grid->fieldPtr_[ID]; 

    // store incremented value at each point 
    int i,j,k; 
    int iter=-1; 
    for (i = 0; i<grid->nxTot_; ++i) { 
        for (j=0; j<grid->nyTot_; ++j) { 
            for (k=0; k<grid->nzTot_; ++k) { 
                field[i][j][k] = ++iter; 
            } 
        } 
    } 

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
    int ydir=1; 
    int zdir=2; 

    int dexYLeftGhost = 0; 
    int dexYRightGhost = grid->fieldSize_[type][ydir] - 1; 
    int dexYLeftReal = ++dexYLeftGhost; 
    int dexYRightReal = --dexYRightGhost; 

    int dexZLeftGhost = 0; 
    int dexZRightGhost = grid->fieldSize_[type][zdir] - 1; 
    int dexZLeftReal = ++dexZLeftGhost; 
    int dexZRightReal = --dexZRightGhost; 

    int val = 0; 
    for (i = 0; i<grid->nxTot_; ++i) { 
        for (j=0; j<grid->nyTot_; ++j) { 
            for (k=0; k<grid->nzTot_; ++k) { 
                val = field[i][j][k]; 
                
                // sums for constant y 
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

                // sums for constant z 
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

    for (i = 0; i<grid->nxTot_; ++i) { 
        for (j=0; j<grid->nyTot_; ++j) { 
            for (k=0; k<grid->nzTot_; ++k) { 
                val = field[i][j][k]; 
                
                // newSums for constant y 
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

                // newSums for constant z 
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

    // check that BCs were updated appropriately 
    EXPECT_EQ(newSumYLeftGhost,sumYRightReal); 
    EXPECT_EQ(newSumYRightGhost,sumYLeftReal); 
    EXPECT_EQ(newSumZLeftGhost,sumZRightReal); 
    EXPECT_EQ(newSumZRightGhost,sumZLeftReal); 

    // check that real cells were unaffected 
    EXPECT_EQ(newSumYLeftReal,sumYLeftReal); 
    EXPECT_EQ(newSumYRightReal,sumYRightReal); 
    EXPECT_EQ(newSumZLeftReal,sumZLeftReal); 
    EXPECT_EQ(newSumYRightReal,sumYRightReal); 
}




