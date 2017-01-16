#define MAIN_CPP
#include <stdio.h>
#include "gtest/gtest.h"
#include "math.h"
#include "poisson.hpp"
#include "particle_handler.hpp"
#include "particle_bc_factory.hpp"
#include "field_bc_factory.hpp"

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

// tests requiring internal examination
class DepositJTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    //MPI
    rank_MPI = 0;
    size_MPI = 1;
    //Read input file
    Input *input =  new Input();
    char filename[200];
    sprintf(filename, "test_data/depositJ_unittest_data.txt");
    input->readinfo(filename);
    input_info = input->getinfo();

    //Initialize the domain, particle handler, boundary conditions and grid.
    domain = new Domain(input_info);
    part_handler = new Particle_Handler();
    //part_handler->setPusher(new Boris());
    bc = Part_BC_Factory::getInstance().constructConditions(domain,input_info);
    //part_handler->setParticleBoundaries(bc);
    grid = new Poisson_Solver(domain, input_info);
    Field_BC_Factory::getInstance().Construct(domain,grid,input_info);
    part_handler->Load(input_info,domain);    
  }

  virtual void TearDown() {
    delete grid;
  }

  Poisson_Solver *grid;
  Particle_Handler *part_handler;
  Domain *domain;
  Input_Info_t *input_info;
  BC_Particle** bc;
};

// Test particle current-deposition is working.
TEST_F(DepositJTest, sumOverJandRho) {
  int np = part_handler->nParticles();
  double mod_qv;
  double q1;
  double q2;

  if (np != 2) {
    printf("This test did not fail, but can only be run for 2 particles. Change input file and rerun.");
  } else {
    //Get full vector of particles.
    std::vector<Particle> parts = *part_handler->getParticleVector();

    //Calculate particle |v| and q
    q1 = parts[0].q;
    q2 = parts[1].q;
    mod_qv = pow(
		 pow(q1*parts[0].v[0]+q2*parts[1].v[0],2.0)+
		 pow(q1*parts[0].v[1]+q2*parts[1].v[1],2.0)+
		 pow(q1*parts[0].v[2]+q2*parts[1].v[2],2.0)
	     ,0.5);

    //Perform the deposition of the current and charge (note: depositRho=true)
    part_handler->depositRhoJ(grid,true,domain,input_info);

    //Get cell volume
    double lcell[3] = {};
    for (int i=0; i<3; i++) lcell[i] = grid->getStepSize(i);
    double cell_volume = lcell[0]*lcell[1]*lcell[2];

    // sum up all cells' J
    int i,j,k;
    double JxSum = 0.0;
    double JySum = 0.0;
    double JzSum = 0.0;
    double rhoSum = 0.0;
    for (i=0; i<grid->nxTot_; ++i) {
      for (j=0; j<grid->nyTot_; ++j) {
	for (k=0; k<grid->nzTot_; ++k) {
	  JxSum += grid->Jx_[i][j][k];
	  JySum += grid->Jy_[i][j][k];
	  JzSum += grid->Jz_[i][j][k];
	  rhoSum += grid->rho_[i][j][k];
	}
      }
    }

    //Calculate deposited |J|, and expected calculated |J|
    double modJ = pow(pow(JxSum,2.0)+pow(JySum,2.0)+pow(JzSum,2.0),0.5);
    double calcModJ = UNIT_RHOJ * mod_qv/cell_volume;

    //set test tolerance
    double tol = 0.000000000000001;

    //Run the test for |J|
    EXPECT_NEAR(modJ,calcModJ,tol);

    //Calculate deposited |rho|, and expected calculated |rho|
    double mod_rho = rhoSum;
    double calc_mod_rho = UNIT_RHOJ * (q1/cell_volume + q2/cell_volume);

    //Run the test for |rho|
    EXPECT_NEAR(mod_rho,calc_mod_rho,tol);
  }
  
}

// Test field interpolation is working.
TEST_F(DepositJTest, testFieldInterpolation) {
    //Get full vector of particles.
    std::vector<Particle> parts = *part_handler->getParticleVector();

    //Make const E field (Ex,Ey,Ez)=(1.0,2.0,3.0)
    grid->constE(1.0,2.0,3.0);

    //Get particle position vector
    double part_pos[3] = {};
    for (int i=0; i<3; i++) part_pos[i] = parts[0].x[i];

    //Get cell lengths
    double lcell[3] = {};
    for (int i=0; i<3; i++) lcell[i] = grid->getStepSize(i);

    //Get cell field variables: cellvars
    double cellvars[21];
    int pCell = grid->getCellID(part_pos[0],part_pos[1],part_pos[2]);
    grid->getFieldInterpolatorVec(pCell, cellvars);

    //Interpolate fields
    Interpolator *interpolator = new Interpolator();
    interpolator->interpolate_fields(part_pos, lcell, cellvars, &(parts[0].field));

    //Check interpolated field values.
    Field_part *field = new Field_part();
    field = &(parts[0].field);
    EXPECT_NEAR(field->e1,1.0,0.00001);
    EXPECT_NEAR(field->e2,2.0,0.00001);
    EXPECT_NEAR(field->e3,3.0,0.00001);

}

//Test Poisson solver is working.
TEST_F(DepositJTest, testPoissonSolver) {
  // Deposit charge and current from particles to grid
  //part_handler->depositRhoJ(grid,true,domain,input_info);
  grid->constJ(0.0,0.0,0.0);
  grid->constRho(0.0);
  grid->rho_[2][2][2] = 1.0;
  grid->rho_[4][4][4] = -1.0;

  // sum charge and current on MPI boundaries
  grid->executeBC(-2); // boundary for Rho,J

  // initialize fields (POISSON SOLVE)
  grid->InitializeFields(input_info);

  ////Define constants used to iterate Poisson's equation
  double dx_ = grid->dx_;
  double dy_ = grid->dy_;
  double dz_ = grid->dz_;
  double celldist2 = pow(dx_, 2.0)*pow(dy_, 2.0) + pow(dy_, 2.0)*pow(dz_, 2.0) + pow(dx_, 2.0)*pow(dz_, 2.0);
  double ax = pow(dy_, 2.0) * pow(dz_, 2.0) / (2.0 * celldist2);
  double ay = pow(dx_, 2.0) * pow(dz_, 2.0) / (2.0 * celldist2);
  double az = pow(dx_, 2.0) * pow(dy_, 2.0) / (2.0 * celldist2);
  double af = pow(dx_, 2.0) * pow(dy_, 2.0) * pow(dz_, 2.0) / (2.0 * celldist2);

  // limits
  double iBeg_ = grid->iBeg_;
  double jBeg_ = grid->jBeg_;
  double kBeg_ = grid->kBeg_;
  double nxReal_ = grid->nxReal_;
  double nyReal_ = grid->nyReal_;
  double nzReal_ = grid->nzReal_;
  int iEnd = nxReal_ + iBeg_;
  int jEnd = nyReal_ + jBeg_;
  int kEnd = nzReal_ + kBeg_;

  double poisson_step = 0.0;
  double*** u1 = grid->phi1_;
  for ( int i=iBeg_; i<iEnd; i++ ) {
    for ( int j=jBeg_; j<jEnd; j++ ) {
      for ( int k=kBeg_; k<kEnd; k++ ) {
	//Calculate poisson step
	poisson_step = u1[i][j][k] - 
	             ( ax*(u1[i-1][j][k]+u1[i+1][j][k])
	             + ay*(u1[i][j-1][k]+u1[i][j+1][k])
	             + az*(u1[i][j][k-1]+u1[i][j][k+1])
		     - af*grid->rho_[i][j][k]*(-4*M_PI));

	EXPECT_NEAR(0.0,poisson_step,0.000000001);
      }
    }
  }
}
