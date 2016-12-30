#include "gtest/gtest.h"
#include "input.hpp"

TEST(ReadinputTest, ValuesCorrect) {
  Input_Info_t input_info;
  
  char filename[100];

  sprintf(filename, "data/test.txt");

  readinput(filename, &input_info, 1);

  EXPECT_EQ(8, input_info.nCell[0]);
  EXPECT_EQ(8, input_info.nCell[1]);
  EXPECT_EQ(8, input_info.nCell[2]);

  EXPECT_EQ(16, input_info.np);
}


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
