#include "hdf5io.hpp"
#include <assert.h>
#include <stdlib.h>
#include <string.h>

Hdf5IO::Hdf5IO(const char* filename, Grid* grid, Domain* domain, const int which_fields) 
	: which_fields_(which_fields),
  	  nProcxyz_(domain->getnProcxyz()),
	  myijk_(domain->getmyijk())
{
  // set up property lists
  // file access
  file_access_plist_ = H5Pcreate(H5P_FILE_ACCESS);
#if USE_MPI
  H5Pset_fapl_mpio(file_access_plist_, MPI_COMM_WORLD, MPI_INFO_NULL);
#else
  file_access_plist_ = H5P_DEFAULT;
#endif
  // data transfer
  data_xfer_plist_ = H5Pcreate(H5P_DATASET_XFER);
#if USE_MPI
  H5Pset_dxpl_mpio(data_xfer_plist_, H5FD_MPIO_COLLECTIVE);
#else
  data_xfer_plist_ = H5P_DEFAULT;
#endif
  assert(file_access_plist_>=0);
  assert(data_xfer_plist_>=0);

  // open file
  file_id_ = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, file_access_plist_);
  assert(file_id_>=0);

  // initialize fields diagnostics
  if(which_fields_>=0) {
    // create field group in hdf5 file
    fields_group_id_ = H5Gcreate2(file_id_, "/fields", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    assert(fields_group_id_>=0);
  }

  int* nxyzTot = new int[3];
  grid->getnxyzTot(nxyzTot);
  int* phys_dims = new int[3];

  if(which_fields_==ALL || which_fields_==E) {
    grid->getDimPhys(grid->getExID(), phys_dims);
    Ex_tsio_ = new FieldTimeseriesIO(this, phys_dims, nxyzTot, "Ex");

    grid->getDimPhys(grid->getEyID(), phys_dims);
    Ey_tsio_ = new FieldTimeseriesIO(this, phys_dims, nxyzTot, "Ey");

    grid->getDimPhys(grid->getEzID(), phys_dims);
    Ez_tsio_ = new FieldTimeseriesIO(this, phys_dims, nxyzTot, "Ez");
  }
  if(which_fields_==ALL || which_fields_==B) {
    grid->getDimPhys(grid->getBxID(), phys_dims);
    Bx_tsio_ = new FieldTimeseriesIO(this, phys_dims, nxyzTot, "Bx");

    grid->getDimPhys(grid->getByID(), phys_dims);
    By_tsio_ = new FieldTimeseriesIO(this, phys_dims, nxyzTot, "By");

    grid->getDimPhys(grid->getBzID(), phys_dims);
    Bz_tsio_ = new FieldTimeseriesIO(this, phys_dims, nxyzTot, "Bz");
  }
  if(which_fields_==ALL || which_fields_==J) {
    grid->getDimPhys(grid->getJxID(), phys_dims);
    Jx_tsio_ = new FieldTimeseriesIO(this, phys_dims, nxyzTot, "Jx");

    grid->getDimPhys(grid->getJyID(), phys_dims);
    Jy_tsio_ = new FieldTimeseriesIO(this, phys_dims, nxyzTot, "Jy");

    grid->getDimPhys(grid->getJzID(), phys_dims);
    Jz_tsio_ = new FieldTimeseriesIO(this, phys_dims, nxyzTot, "Jz");
  }
  if(which_fields_==ALL || which_fields_==rho) {
    grid->getDimPhys(grid->getrhoID(), phys_dims);
    rho_tsio_ = new FieldTimeseriesIO(this, phys_dims, nxyzTot, "rho");
  }
}

Hdf5IO::~Hdf5IO() {
  // close property lists
  H5Pclose(file_access_plist_);
  H5Pclose(data_xfer_plist_);

  // close fields group
  H5Gclose(fields_group_id_);

  // close file
  H5Fclose(file_id_);
}

int Hdf5IO::writeFields(Grid* grid) {
  double**** fieldPtr = grid->getFieldPtr();  
  assert(fieldPtr!=NULL);
  if(which_fields_==ALL || which_fields_==E) {
    Ex_tsio_->writeField(fieldPtr[grid->getExID()]);
    Ey_tsio_->writeField(fieldPtr[grid->getEyID()]);
    Ez_tsio_->writeField(fieldPtr[grid->getEzID()]);
  }
  if(which_fields_==ALL || which_fields_==B) {
    Bx_tsio_->writeField(fieldPtr[grid->getBxID()]);
    By_tsio_->writeField(fieldPtr[grid->getByID()]);
    Bz_tsio_->writeField(fieldPtr[grid->getBzID()]);
  }
  if(which_fields_==ALL || which_fields_==J) {
    Jx_tsio_->writeField(fieldPtr[grid->getJxID()]);
    Jy_tsio_->writeField(fieldPtr[grid->getJyID()]);
    Jz_tsio_->writeField(fieldPtr[grid->getJzID()]);
  }
  if(which_fields_==ALL || which_fields_==rho) {
    rho_tsio_->writeField(fieldPtr[grid->getrhoID()]);
  }

  return 0;
}

FieldTimeseriesIO::FieldTimeseriesIO(Hdf5IO* io, int* phys_dims, int* nxyzTot, std::string fieldname) 
	: ndims_(4),
  	  nProcxyz_(io->getnProcxyz()),
	  myijk_(io->getmyijk()),
	  file_id_(io->getFileID()),
	  data_xfer_plist_(io->getDataXferPlist()),
	  fields_group_id_(io->getFieldsGroupID())
{
  
  // allocate some (temporary) arrays needed in the loop below
  hsize_t* filespace_dims = new hsize_t[ndims_];
  hsize_t* filespace_maxdims = new hsize_t[ndims_];
  hsize_t* memspace_dims = new hsize_t[ndims_];

  // allocating array storing block dimensions
  field_block_ = new hsize_t[ndims_]; 
  assert(field_block_!=NULL);

  // use physical dimensions for block spatial dimensions
  field_block_[0] = phys_dims[0];
  field_block_[1] = phys_dims[1];
  field_block_[2] = phys_dims[2];
  // each block is a snapshot in time, so set time dimension to 1
  field_block_[3] = 1;

  // dimensions of 4D field data in file
  // each proc will write a field block, so need field_block*nProc in each dimension
  filespace_dims[0] = filespace_maxdims[0] = field_block_[0]*nProcxyz_[0];
  filespace_dims[1] = filespace_maxdims[1] = field_block_[1]*nProcxyz_[1];
  filespace_dims[2] = filespace_maxdims[2] = field_block_[2]*nProcxyz_[2];
  // initialize time dimension to 0. we will extend before each time write
  filespace_dims[3] = 0;
  // let time dimension be extensible
  filespace_maxdims[3] = H5S_UNLIMITED;

  // create a dataspace for the file
  filespace_ = H5Screate_simple(ndims_, filespace_dims, filespace_maxdims);
  assert(filespace_>=0);

  // set up file dataspace to be chunked
  hid_t create_chunks_plist = H5Pcreate(H5P_DATASET_CREATE);
  // make chunks the size of the field blocks
  H5Pset_chunk(create_chunks_plist, ndims_, field_block_);
  assert(create_chunks_plist>=0);

  // create a dataset for each field
  field_dataset_ = H5Dcreate(fields_group_id_, fieldname.c_str(), H5T_NATIVE_DOUBLE, 
                             filespace_, H5P_DEFAULT, create_chunks_plist, H5P_DEFAULT);
  assert(field_dataset_>=0);
  H5Pclose(create_chunks_plist); 

  // above we have prepared the hdf5 file for writing.
  // now we need to define the dataspace of each field in the program memory
  memspace_dims[0] = nxyzTot[0];
  memspace_dims[1] = nxyzTot[1];
  memspace_dims[2] = nxyzTot[2];
  memspace_dims[3] = 1;
  memspace_ = H5Screate_simple(ndims_, memspace_dims, NULL);
  assert(memspace_>=0);
  
  delete [] filespace_dims;
  delete [] filespace_maxdims;
  delete [] memspace_dims;
}

FieldTimeseriesIO::~FieldTimeseriesIO() {
  H5Sclose(memspace_);
  H5Sclose(filespace_);
  H5Dclose(field_dataset_);
  delete[] field_block_;
}

int FieldTimeseriesIO::writeField(double*** field_data) {
  int status = -1;

  // Note: field_block_, filespace_, memspace_ allocated and set in constructor

  hsize_t* offset = new hsize_t[ndims_];
  hsize_t* stride = new hsize_t[ndims_];
  hsize_t* count = new hsize_t[ndims_];
  assert(offset!=NULL);
  assert(stride!=NULL);
  assert(count!=NULL);

  // extend filespace and dataset for this iwrite
  hsize_t* currdims = new hsize_t[ndims_];
  hsize_t* maxdims = new hsize_t[ndims_];
  hsize_t* newdims = new hsize_t[ndims_];
  H5Sget_simple_extent_dims(filespace_, currdims, maxdims);
  newdims[0] = currdims[0];
  newdims[1] = currdims[1];
  newdims[2] = currdims[2];
  newdims[3] = currdims[3]+1;
  status = H5Sset_extent_simple(filespace_, ndims_, newdims, maxdims);
  status = H5Dset_extent(field_dataset_, newdims);
  assert(status>=0);

  // select the subset of the file dataspace that this proc will be writing to
  offset[0] = myijk_[0]*field_block_[0];
  offset[1] = myijk_[1]*field_block_[1];
  offset[2] = myijk_[2]*field_block_[2];
  offset[3] = newdims[3]-1;
  stride[0] = stride[1] = stride[2] = stride[3] = 1;
  count[0] = count[1] = count[2] = count[3] = 1;
  
  status = H5Sselect_hyperslab(filespace_, H5S_SELECT_SET, 
				   offset, stride, count, field_block_);
  assert(status>=0);

  // select the subset of the memory dataspace that we are writing from
  // skip first cell in each spatial dim, which is ghost
  offset[0] = 1;
  offset[1] = 1;
  offset[2] = 1;
  offset[3] = 0;
  stride[0] = stride[1] = stride[2] = stride[3] = 1;
  count[0] = count[1] = count[2] = count[3] = 1;

  status = H5Sselect_hyperslab(memspace_, H5S_SELECT_SET, 
				   offset, stride, count, field_block_);
  assert(status>=0);

  // write from field_data in memspace to filespace in file
  status = H5Dwrite(field_dataset_, H5T_NATIVE_DOUBLE, 
		      memspace_, filespace_,
		      data_xfer_plist_, field_data[0][0]); 
  assert(status>=0);

  delete [] offset;
  delete [] stride;
  delete [] count;
  
  return status;
}



