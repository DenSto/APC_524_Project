#include "gtest/gtest.h"
#include "grid.hpp"

TEST(GridTestContiguousAllocation, doesNotCrash) {
  int nxyz [3] = {3,4,5};
  double xyz0 [3] = {0,0,0};
  double Lxyz [3] = {1,1,1};
  int nGhosts = 1; 

  Grid *grid = new Grid(nxyz, nGhosts, xyz0, Lxyz);
  EXPECT_EQ(60, grid->getNumberOfCells());
  delete grid;
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
