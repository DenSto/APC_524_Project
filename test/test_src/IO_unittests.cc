#include "gtest/gtest.h"
#include "IO.hpp"

TEST(ReadinputTest, ValuesCorrect) {
  Input_Info_t input_info;
  
  char filename[100];

  sprintf(filename, "dummy");

  readinput(filename, &input_info);

  EXPECT_EQ(4, input_info.nx);
  EXPECT_EQ(8, input_info.np);
}


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
