#include "gtest/gtest.h"
#include "grid.hpp"

TEST(GridGetNumberOfCellsTest, returnsCorrectNumberOfCells) {
  int nxyz [3] = {3,4,5};
  double xyz0 [3] = {0,0,0};
  double Lxyz [3] = {1,1,1};

  Grid *grid = new Grid(nxyz, 1, xyz0, Lxyz);
  EXPECT_EQ(60, grid->getNumberOfCells());
  delete grid;
}

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
