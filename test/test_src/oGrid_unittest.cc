#include "gtest/gtest.h"
#include "grid.hpp"

TEST(GridGetNumberOfCellsTest, returnsCorrectNumberOfCells) {
  int nxyz [3] = {3,4,5};
  double xyz0 [3] = {0,0,0};
  double Lxyz [3] = {1,1,1};

  Grid *grid = new Grid(nxyz, xyz0, Lxyz);
  EXPECT_EQ(60, grid->getNumberOfCells());
}

class OGridTest : public ::testing::Test {
   protected:
      virtual void SetUp() {
         int nxyz [3] = {5,10,20};
         double xyz0 [3] = {0,0,0};
         double Lxyz [3] = {1,1,1};

         grid = new Grid(nxyz, xyz0, Lxyz);

      }

      virtual void TearDown() {
         delete grid;
      }

      Grid *grid;

};

TEST_F(OGridTest, getCellIDReturnsInternalCells) {
  double xyz [3] = { .501, .501, .501};

  EXPECT_EQ(510, grid->getCellID(xyz));
}

TEST_F(OGridTest, getCellIDReturnsInternalCells) {
  double xyz [3] = { .501, .501, .501};

  EXPECT_EQ(510, grid->getCellID(xyz));
}




int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
