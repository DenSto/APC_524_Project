#include "gtest/gtest.h"
#include "grid.hpp"

TEST(GridGetNumberOfCellsTest, returnsCorrectNumberOfCells) {
  int nxyz [3] = {3,4,5};
  double xyz0 [3] = {0,0,0};
  double Lxyz [3] = {1,1,1};

  Grid *grid = new Grid(nxyz, xyz0, Lxyz);
  EXPECT_EQ(60, grid->getNumberOfCells());
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
