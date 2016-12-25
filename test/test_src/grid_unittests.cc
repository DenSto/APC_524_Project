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
         int nxyz [3] = {3,5,7};
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

// test that fieldSize is correct  
TEST_F(gridPrivateTest, fieldSizeTest) {
    int xdir=0; 
    int ydir=1; 
    int zdir=2; 

    // test a set of edges
    EXPECT_EQ(grid->fieldSize_[grid->edgeXID_][xdir],grid->nx_); 
    EXPECT_EQ(grid->fieldSize_[grid->edgeXID_][ydir],grid->nyTot_); 
    EXPECT_EQ(grid->fieldSize_[grid->edgeXID_][zdir],grid->nzTot_); 

    // test a set of faces 
    EXPECT_EQ(grid->fieldSize_[grid->faceYID_][xdir],grid->nx_); 
    EXPECT_EQ(grid->fieldSize_[grid->faceYID_][ydir],grid->nyTot_); 
    EXPECT_EQ(grid->fieldSize_[grid->faceYID_][zdir],grid->nz_);
    
    // test a set of vertices
    EXPECT_EQ(grid->fieldSize_[grid->vertID_][xdir],grid->nxTot_); 
    EXPECT_EQ(grid->fieldSize_[grid->vertID_][ydir],grid->nyTot_); 
    EXPECT_EQ(grid->fieldSize_[grid->vertID_][zdir],grid->nzTot_); 
}

// test that fieldPtr works 
TEST_F(gridPrivateTest, fieldPtrTest) {
    double Exval = 1; 
    double Byval = 2; 
    double rhoval = 3; 

    int i,j,k; 
    for (i=0; i<grid->nxTot_; ++i) { 
        for (j=0; j<grid->nyTot_; ++j) { 
            for (k=0; k<grid->nzTot_; ++k) { 
                grid->Ex_[i][j][k] = Exval; 
                grid->By_[i][j][k] = Byval; 
                grid->rho_[i][j][k] = rhoval; 
            }
        }
    }

    double*** ExPtr = grid->fieldPtr_[grid->ExID_]; 
    double*** ByPtr = grid->fieldPtr_[grid->ByID_];
    double*** rhoPtr = grid->fieldPtr_[grid->rhoID_];

    EXPECT_EQ(grid->Ex_[1][1][1],ExPtr[1][1][1]); 
    EXPECT_EQ(grid->By_[1][1][1],ByPtr[1][1][1]); 
    EXPECT_EQ(grid->rho_[1][1][1],rhoPtr[1][1][1]); 
}


// test that zeroField methods work 
TEST_F(gridPrivateTest,zeroFields) {
 
    // zero all fields 
    grid->zeroE(); 
    grid->zeroB(); 
    grid->zeroJ(); 
    grid->zeroRho(); 

    // sum each field 
    double ExSum = sumField(grid->Ex_); 
    double BySum = sumField(grid->By_); 
    double JzSum = sumField(grid->Jz_); 
    double rhoSum = sumField(grid->rho_); 

    // verify that each field sums to 0
    EXPECT_EQ(ExSum,0); 
    EXPECT_EQ(BySum,0); 
    EXPECT_EQ(JzSum,0); 
    EXPECT_EQ(rhoSum,0); 

    // set one point in Ex to this value 
    double addval = -2.0; 
    grid->Ex_[1][1][1] += addval; 
    ExSum = sumField(grid->Ex_); 
    // verify that the sum is no longer zero 
    EXPECT_EQ(ExSum,addval); 
    // zero E again 
    grid->zeroE(); 
    // verify that E again sums to zero 
    ExSum = sumField(grid->Ex_); 
    EXPECT_EQ(ExSum,0); 

}




