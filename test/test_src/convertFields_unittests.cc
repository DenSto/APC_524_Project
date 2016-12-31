#include "gtest/gtest.h"
#include "poisson.hpp"
#include "math.h"

// Tests of input/output to public methods
class ConvertTest : public ::testing::Test {
   protected:
      virtual void SetUp() {
         int nxyz [3] = {3,5,7};
         int nGhosts = 1; 
         double xyz0 [3] = {0,0,0};
         double Lxyz [3] = {1,1,1};

         grid = new Poisson_Solver(nxyz, nGhosts, xyz0, Lxyz);

      }

      virtual void TearDown() {
         delete grid;
      }

      Poisson_Solver *grid;

};

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}


// tests requiring internal examination
class ConvertPrivateTest : public ::testing::Test {
   protected:
      virtual void SetUp() {
         int nxyz [3] = {3,5,7};
         int nGhosts = 1; 
         double xyz0 [3] = {0,0,0};
         double Lxyz [3] = {1,1,1};

         grid = new Poisson_Solver(nxyz, nGhosts, xyz0, Lxyz);

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

      Poisson_Solver *grid;

};

// test that conversion phi to E is correct
// case: constant phi 
TEST_F(ConvertPrivateTest, constantPhiTest) {

    // set phi1 to constant 
    int i,j,k; 
    double phival = 1.0; 
    for (i=0; i<grid->nxTot_; ++i) { 
        for (j=0; j<grid->nyTot_; ++j) { 
            for (k=0; k<grid->nzTot_; ++k) { 
                grid->phi1_[i][j][k] = phival; 
            } 
        } 
    } 

    // derive E 
    grid->zeroE(); 
    grid->phiToE(); 

    // for constant phi, E should vanish everywhere 
    double  ExSum = sumField(grid->Ex_); 
    EXPECT_EQ(ExSum,0); 
}

