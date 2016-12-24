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

TEST_F(GridTest, getStepSize) {
    // placeholder test
    EXPECT_EQ(1,1);
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}


// tests requiring internal examination
class gridPrivateTest : public ::testing::Test {
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

// test that the fieldSizes are correct 
TEST_F(gridPrivateTest, fieldSizeTest) {
    int xdir=0; 
    int ydir=1; 
    int zdir=2; 

    EXPECT_EQ(grid->fieldSize_[grid->edgeXID_][xdir],grid->nx_); 
    EXPECT_EQ(grid->fieldSize_[grid->edgeXID_][ydir],grid->nyTot_); 
    EXPECT_EQ(grid->fieldSize_[grid->edgeXID_][zdir],grid->nzTot_); 

    EXPECT_EQ(grid->fieldSize_[grid->edgeYID_][xdir],grid->nxTot_); 
    EXPECT_EQ(grid->fieldSize_[grid->edgeYID_][ydir],grid->ny_); 
    EXPECT_EQ(grid->fieldSize_[grid->edgeYID_][zdir],grid->nzTot_);

    EXPECT_EQ(grid->fieldSize_[grid->edgeZID_][xdir],grid->nxTot_); 
    EXPECT_EQ(grid->fieldSize_[grid->edgeZID_][ydir],grid->nyTot_); 
    EXPECT_EQ(grid->fieldSize_[grid->edgeZID_][zdir],grid->nz_); 

    EXPECT_EQ(grid->fieldSize_[grid->faceXID_][xdir],grid->nxTot_); 
    EXPECT_EQ(grid->fieldSize_[grid->faceXID_][ydir],grid->ny_); 
    EXPECT_EQ(grid->fieldSize_[grid->faceXID_][zdir],grid->nz_); 

    EXPECT_EQ(grid->fieldSize_[grid->faceYID_][xdir],grid->nx_); 
    EXPECT_EQ(grid->fieldSize_[grid->faceYID_][ydir],grid->nyTot_); 
    EXPECT_EQ(grid->fieldSize_[grid->faceYID_][zdir],grid->nz_);

    EXPECT_EQ(grid->fieldSize_[grid->faceZID_][xdir],grid->nx_); 
    EXPECT_EQ(grid->fieldSize_[grid->faceZID_][ydir],grid->ny_); 
    EXPECT_EQ(grid->fieldSize_[grid->faceZID_][zdir],grid->nzTot_); 

    EXPECT_EQ(grid->fieldSize_[grid->vertID_][xdir],grid->nxTot_); 
    EXPECT_EQ(grid->fieldSize_[grid->vertID_][ydir],grid->nyTot_); 
    EXPECT_EQ(grid->fieldSize_[grid->vertID_][zdir],grid->nzTot_); 
}

// TO DO: test that the fieldPtrs actually point to their fields
// (test: fill fields with different values, test that direct access 
// and access via pointer agree)



