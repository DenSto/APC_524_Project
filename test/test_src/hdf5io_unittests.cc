#define MAIN_CPP
#include "gtest/gtest.h"
#include "domain.hpp"
#include "grid.hpp"
#include "hdf5io.hpp"
#include "math.h"

// Tests of input/output to public methods
class FieldIOTest : public ::testing::Test {
protected:
  virtual void SetUp() {

  }

  virtual void TearDown() {
    delete grid;
  }

  Domain *domain;
  Grid *grid;
  Hdf5IO *hdf5io;
};

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

TEST_F(FieldIOTest, writeField) {
    int nxyz [3] = {3,4,5};
    int nProc [3] = {1,1,1};
    int nGhosts = 1; 
    double xyz0 [3] = {0,0,0};
    double Lxyz [3] = {1,1,1};

    int which_fields = 1; //only write E
    int nwrite = 0;

    size_MPI = 1;
    rank_MPI = 0;
    domain = new Domain(nxyz, nProc, xyz0, Lxyz);
    grid = new Grid(nxyz, nGhosts, xyz0, Lxyz);
    hdf5io = new Hdf5IO("test_data/hdf5io_unittests.h5", grid, domain, which_fields);

    int i,j,k; 
    for (i=0; i<grid->nx_+1; ++i) { 
        for (j=0; j<grid->ny_+1; ++j) { 
            for (k=0; k<grid->nz_+1; ++k) { 
                grid->Ex_[i][j][k] = -( k + 10*j + 100*i ); 
                grid->Ey_[i][j][k] = -10*( k + 10*j + 100*i ); 
                grid->Ez_[i][j][k] = -100*( k + 10*j + 100*i ); 
            }
        }
    }
    hdf5io->writeFields(grid);
    nwrite++;

    // change Ex
    int dim_phys[3];
    grid->getDimPhys(grid->getExID(), dim_phys); 

    int iBeg = 1;
    int jBeg = 1;
    int kBeg = 1;
    int iEnd = dim_phys[0]+1;
    int jEnd = dim_phys[1]+1;
    int kEnd = dim_phys[2]+1;

    for (i=iBeg; i<iEnd; ++i) { 
        for (j=jBeg; j<jEnd; ++j) { 
            for (k=kBeg; k<kEnd; ++k) { 
                grid->Ex_[i][j][k] = (k) + 10*(j) + 100*(i); 
            }
        }
    }
    hdf5io->writeFields(grid);
    nwrite++;

    // check file
    hid_t file_id = hdf5io->getFileID();
    ASSERT_GE(file_id, 0);
    hid_t dset_id;
    dset_id = H5Dopen2(file_id, "/fields/Ey", H5P_DEFAULT);
    ASSERT_GE(dset_id, 0);
    dset_id = H5Dopen2(file_id, "/fields/Ez", H5P_DEFAULT);
    ASSERT_GE(dset_id, 0);
    // check Ex values
    dset_id = H5Dopen2(file_id, "/fields/Ex", H5P_DEFAULT);
    ASSERT_GE(dset_id, 0);
    hid_t filespace = H5Dget_space(dset_id);
    ASSERT_GE(filespace, 0);

    hsize_t* dims = new hsize_t[4];
    hsize_t* maxdims = new hsize_t[4];
    H5Sget_simple_extent_dims(filespace, dims, maxdims);

    ASSERT_EQ(dims[0], dim_phys[0]);
    ASSERT_EQ(dims[1], dim_phys[1]);
    ASSERT_EQ(dims[2], dim_phys[2]);
    ASSERT_EQ(dims[3], nwrite);

    double* check_data = new double[dim_phys[0]*dim_phys[1]*dim_phys[2]*nwrite];
    H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, filespace, H5P_DEFAULT, check_data);

    for (i=iBeg; i<iEnd; ++i) { 
        for (j=jBeg; j<jEnd; ++j) { 
            for (k=kBeg; k<kEnd; ++k) { 
                   EXPECT_EQ(-grid->Ex_[i][j][k], check_data[0 + nwrite*(k-1) + nwrite*dim_phys[2]*(j-1) + nwrite*dim_phys[2]*dim_phys[1]*(i-1)]);
                   EXPECT_EQ(grid->Ex_[i][j][k], check_data[1 + nwrite*(k-1) + nwrite*dim_phys[2]*(j-1) + nwrite*dim_phys[2]*dim_phys[1]*(i-1)]);
            }
        }
    }
    

}

