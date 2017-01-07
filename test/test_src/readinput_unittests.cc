#define MAIN_CPP
#include "gtest/gtest.h"
#include "input.hpp"

TEST(ReadinputTest, ValuesCorrect) {
  Input *input =  new Input();
  char filename[100];
  sprintf(filename, "test_data/readinput_unittests.txt");
  input->readinfo(filename);
  Input_Info_t *input_info = input->getinfo();

  EXPECT_EQ(8, input_info->nCell[0]);
  EXPECT_EQ(8, input_info->nCell[1]);
  EXPECT_EQ(8, input_info->nCell[2]);

  EXPECT_EQ(16, input_info->np);
}


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
