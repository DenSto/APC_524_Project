#define MAIN_CPP
#include "gtest/gtest.h"
#include "particle.hpp"
#include "particle_handler.hpp"
#include "boris.hpp"
#include "relativisticBoris.hpp"
#include <time.h>
#include "RNG.hpp"
#include "math.h"

// Tests of input/output to public methods
class ParticleTest : public ::testing::Test {
protected:
  virtual void SetUp() {
	srand (time(NULL));
	RNG = new Random_Number_Generator(rand());
	dt = 0.01;
  }

  virtual void TearDown() {
	  delete RNG;
  }

  Random_Number_Generator* RNG;
  double dt;
};

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

// test that a particle can execute Larmor motion
TEST_F(ParticleTest,larmorTest){
	Relativistic_Boris* boris = new Relativistic_Boris(); 
	Particle p = new_particle();
	Field_part f = new_particle_field();

	p.q = 1;
	p.x[0] = 1.0;
	p.v[1] = 1.0;

	f.b3 = 1.0;

	for(int i = 0; i < 100000; i++){
		boris->Step(&p,&f,dt);
	}

	double rsqr= p.x[0]*p.x[0] + p.x[1]*p.x[1] + p.x[2]+p.x[2];
	double r = sqrt(rsqr);

	EXPECT_NEAR(r,1.0,0.1); // Gyroradius of 1 not quite centered at 
							// 0, since velocities are at half-times

	delete boris;
}

// test that a particle cannot travel faster than speed of light
TEST_F(ParticleTest,lightSpeedTest){
	Relativistic_Boris* boris = new Relativistic_Boris(); 
	Particle p = new_particle();

	p.q = 1.0;

	Field_part f = new_particle_field();

	double fMag = 5;

	f.e1=fMag*(2*RNG->getUniform()-1.0);
	f.e2=fMag*(2*RNG->getUniform()-1.0);
	f.e3=fMag*(2*RNG->getUniform()-1.0);

	for(int i = 0; i < 100000; i++){
		boris->Step(&p,&f,dt);
	}
	
	double vsqr= p.v[0]*p.v[0] + p.v[1]*p.v[1] + p.v[2]*p.v[2];

	
	ASSERT_TRUE(1.0-vsqr > 0); // Can't be larger than one...
	EXPECT_NEAR(1.0-sqrt(vsqr),0.0,1e-4); // But should be close after to one after constant acceleration

	delete boris;
} 
