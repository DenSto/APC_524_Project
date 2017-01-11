#define MAIN_CPP
#include "gtest/gtest.h"
#include "RNG.hpp"
#include <cmath>
#include "globals.hpp"

#define SAMPLES 100000

TEST(ReadinputTest, ValuesCorrect) {
  Random_Number_Generator *rng = new Random_Number_Generator(1233); //Note seed=1232 gives an outlier >2*std case.

  double aveUni=0, aveG1=0, aveG2=0, varG1=0,varG2=0;
  double shotNoise = 1.0/sqrt((double)SAMPLES);

  double uni1,g1,g2;
  for(int i =0; i < SAMPLES; i++){
 	uni1=rng->getUniform(); 
  	g1=rng->getStandardNormal();
  	g2=rng->getStandardNormal();
	aveUni+=uni1;
	aveG1+=g1;
	aveG2+=g2;

 	varG1+=g1*g1; 
 	varG2+=g2*g2; 
  }

  aveUni/=((double)SAMPLES);
  aveG1/=((double)SAMPLES);
  aveG2/=((double)SAMPLES);
  varG1/=((double)SAMPLES);
  varG2/=((double)SAMPLES);
  
  EXPECT_NEAR(0.5,aveUni,shotNoise);

  EXPECT_NEAR(0,aveG1,2*shotNoise);
  EXPECT_NEAR(0,aveG2,2*shotNoise);
  EXPECT_NEAR(1.0,varG1,2*shotNoise);
  EXPECT_NEAR(1.0,varG2,2*shotNoise);

}


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
