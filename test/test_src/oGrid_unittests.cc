#define MAIN_CPP
#include "gtest/gtest.h"
#include "grid.hpp"
#include "math.h"

TEST(GridGetNumberOfCellsTest, returnsCorrectNumberOfCells) {
  int nxyz [3] = {3,4,5};
  double xyz0 [3] = {0,0,0};
  double Lxyz [3] = {1,1,1};

  Grid *grid = new Grid(nxyz, 1, xyz0, Lxyz);
  EXPECT_EQ(60, grid->getNumberOfCells());
  delete grid;
}

// Tests of input/output to public methods
class OGridTest : public ::testing::Test {
   protected:
      virtual void SetUp() {
         int nxyz [3] = {5,10,20};
         double xyz0 [3] = {0,0,0};
         double Lxyz [3] = {1,1,1};

         grid = new Grid(nxyz, 1, xyz0, Lxyz);

      }

      virtual void TearDown() {
         delete grid;
      }

      Grid *grid;

};

TEST_F(OGridTest, getStepSize) {
  EXPECT_EQ(.2, grid->getStepSize(0));
  EXPECT_EQ(.1, grid->getStepSize(1));
  EXPECT_EQ(.05, grid->getStepSize(2));
}


TEST_F(OGridTest, getCellIDReturnsInternalCells) {
  EXPECT_EQ(510, grid->getCellID(0.501,0.501,0.501));
}

TEST_F(OGridTest, getCellIDReturnsGhostCells) {
  EXPECT_EQ(-1, grid->getCellID( 0, 0, 0));
  EXPECT_EQ(-2, grid->getCellID( 0, 0, 1));
  EXPECT_EQ(-3, grid->getCellID( 0, 0, .5));
  EXPECT_EQ(-4, grid->getCellID( 0, 1, .5));
  EXPECT_EQ(-5, grid->getCellID( 0, .5, .5));
  EXPECT_EQ(-6, grid->getCellID( 1, .5, .5));
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}


// tests requiring internal examination
class oGridInternalTest : public ::testing::Test {
   protected:
      virtual void SetUp() {
         int nxyz [3] = {20,20,200};
         double xyz0 [3] = {0,0,0};
         double Lxyz [3] = {1,1,1};

         grid = new Grid(nxyz, 1, xyz0, Lxyz);

      }

      virtual void TearDown() {
         delete grid;
      }

      Grid *grid;

};
/*
TEST_F(oGridInternalTest, EMWave) {

   // Initialize EM plane wave
   for (int ix = 0; ix < grid -> nx_; ix++) {
      for (int iy = 0; iy < grid -> ny_; iy++) {
         for (int iz = 0; iz < grid -> nz_; iz++) {
            grid -> Ex_[ix][iy][iz] = cos( M_PI * iz);
            grid -> By_[ix][iy][iz] = cos( M_PI * iz + M_PI/2);
         }
      }
   }
   EXPECT_LT( fabs(1 - grid->Ex_[10][10][500]), 1e-12);
   EXPECT_LT( fabs(0 - grid->By_[10][10][500]), 1e-12);

   double T = 2*grid->getStepSize(2);

   double dt = T / 100;

   grid->constJ(0,0,0);
   for (int i = 0; i < 100; i++) {
      grid->evolveFields( dt);
      grid->updatePeriodicGhostCells();
      printf("Ex: %f   By: %f\n", grid->Ex_[10][10][100], grid->By_[10][10][100]);
   }

   EXPECT_LT( fabs(1 - grid->Ex_[5][5][6]), 1e-12);
   EXPECT_LT( fabs(0 - grid->By_[5][5][6]), 1e-12);


}
*/


TEST_F(oGridInternalTest, EMWaveLong) {

   // Initialize EM plane wave
   double dxL = 10; //dx / lambda; choose 10, 20, 50
   for (int ix = 0; ix < grid -> nx_; ix++) {
      for (int iy = 0; iy < grid -> ny_; iy++) {
         for (int iz = 0; iz < grid -> nz_; iz++) {
            grid -> Ex_[ix][iy][iz] = cos( 2 * M_PI * (((double) iz)/dxL));
            grid -> By_[ix][iy][iz] = cos( 2 * M_PI * (((double) iz)/dxL) + M_PI/dxL);
         }
      }
   }
   EXPECT_LT( fabs(1 - grid->Ex_[10][10][100]), 1e-12);
   EXPECT_LT( fabs(cos( M_PI/dxL ) - grid->By_[10][10][100]), 1e-12);
   
   double T = dxL*grid->getStepSize(2);

   double dt = T / 1000;

   for (int i = 0; i <= 250; i++) {
      grid->evolveFields( dt);
      grid->updatePeriodicGhostCells();
      if ( i %10 == 0) {
         //printf("Ex: %f   By: %f \n", grid->Ex_[10][10][100], grid->By_[10][10][500]);
      }
   }
   EXPECT_LT( fabs(0 - grid->Ex_[10][10][100]), 7e-2);
   EXPECT_LT( fabs(cos( - M_PI/2 + M_PI/dxL) - grid->By_[10][10][100]), 7e-2);


   for (int i = 0; i <= 250; i++) {
      grid->evolveFields( dt);
      grid->updatePeriodicGhostCells();
      if ( i %10 == 0) {
         //printf("Ex: %f   By: %f \n", grid->Ex_[10][10][100], grid->By_[10][10][500]);
      }
   }
   EXPECT_LT( fabs(-1 - grid->Ex_[10][10][100]), 7e-2);
   EXPECT_LT( fabs( cos( - M_PI + M_PI/dxL)  - grid->By_[10][10][100]), 7e-2);


   for (int i = 0; i <= 250; i++) {
      grid->evolveFields( dt);
      grid->updatePeriodicGhostCells();
      if ( i %10 == 0) {
         //printf("Ex: %f   By: %f \n", grid->Ex_[10][10][100], grid->By_[10][10][500]);
      }
   }
   EXPECT_LT( fabs(0 - grid->Ex_[10][10][100]), 7e-2);
   EXPECT_LT( fabs(cos( - 3*M_PI/2 + M_PI/dxL)- grid->By_[10][10][100]), 7e-2);


   for (int i = 0; i <= 250; i++) {
      grid->evolveFields( dt);
      grid->updatePeriodicGhostCells();
      if ( i %10 == 0) {
         //printf("Ex: %f   By: %f \n", grid->Ex_[10][10][100], grid->By_[10][10][500]);
      }
   }
   EXPECT_LT( fabs(1 - grid->Ex_[10][10][100]), 7e-2);
   EXPECT_LT( fabs(cos( - 2*M_PI + M_PI/dxL)  - grid->By_[10][10][100]), 7e-2);


}



