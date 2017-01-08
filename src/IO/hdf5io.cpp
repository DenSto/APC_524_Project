#include "hdf5io.hpp"
#include <assert.h>
#include <stdlib.h>
#include <string.h>

Hdf5IO::Hdf5IO(const char* filename) 
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
}

Hdf5IO::~Hdf5IO() {
  // close property lists
  H5Pclose(file_access_plist_);
  H5Pclose(data_xfer_plist_);

  // close file
  H5Fclose(file_id_);
}

FieldTimeseriesIO::FieldTimeseriesIO(Hdf5IO* io, Grid* grid, Domain* domain, const int which_fields, const int totWrites) 
	: nFieldDatasets_(13),
	  totWrites_(totWrites),
  	  ndims_(4),
  	  nProcxyz_(domain->getnProcxyz()),
	  myijk_(domain->getmyijk())
{
  file_id_ = io->getFileID();
  data_xfer_plist_ = io->getDataXferPlist();

  // allocate a dataspace (memory and file) for each field
  memspace_ = (hid_t*) malloc(sizeof(hid_t)*nFieldDatasets_);
  filespace_ = (hid_t*) malloc(sizeof(hid_t)*nFieldDatasets_);
  // allocate a dataset for each field
  field_dataset_ = (hid_t*) malloc(sizeof(hid_t)*nFieldDatasets_);

  // allocate 2D array for storing block dimensions for each field 
  // second dimension allocated in loop below
  field_block_ = new hsize_t*[nFieldDatasets_];
  assert(field_block_!=NULL);
  
  // create field group in hdf5 file
  fields_group_id_ = H5Gcreate2(file_id_, "/fields", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  assert(fields_group_id_>=0);

  char fieldname[8];

  // allocate some (temporary) arrays needed in the loop below
  int* phys_dims = new int[3];
  hsize_t* filespace_dims = new hsize_t[ndims_];
  hsize_t* filespace_maxdims = new hsize_t[ndims_];
  int* nxyzTot = new int[3];
  hsize_t* memspace_dims = new hsize_t[ndims_];

  // this isn't very elegant...
  for(int fieldID=0; fieldID<nFieldDatasets_; fieldID++) {
    //if( !(input->write_E || input->write_all_fields) 
    if( !(which_fields==1 || which_fields==4) 
       && 
      (fieldID==grid->getExID() || fieldID==grid->getEyID() || fieldID==grid->getEzID()) )
    {
      continue;
    }
    //if( !(input->write_J || input->write_all_fields) 
    if( !(which_fields==3 || which_fields==4) 
       && 
      (fieldID==grid->getJxID() || fieldID==grid->getJyID() || fieldID==grid->getJzID()) )
    {
      continue;
    }
    //if( !(input->write_B || input->write_all_fields) 
    if( !(which_fields==2 || which_fields==4) 
       && 
      (fieldID==grid->getBxID() || fieldID==grid->getByID() || fieldID==grid->getBzID() ||
       fieldID==grid->getBx_tm1ID() || fieldID==grid->getBy_tm1ID() || fieldID==grid->getBz_tm1ID()))
    {
      continue;
    }
    //if( !(input->write_rho || input->write_all_fields) 
    if( !(which_fields==0 || which_fields==4) 
       && 
      (fieldID==grid->getrhoID()) )
    {
      continue;
    }

    // finish allocating 2D for storing block dimensions
    field_block_[fieldID] = new hsize_t[ndims_]; 
    assert(field_block_[fieldID]!=NULL);

    // get (local) dimensions of the physical (non-ghost) field region
    grid->getDimPhys(fieldID, phys_dims);

    // use physical dimensions for block spatial dimensions
    field_block_[fieldID][0] = phys_dims[0];
    field_block_[fieldID][1] = phys_dims[1];
    field_block_[fieldID][2] = phys_dims[2];
    // each block is a snapshot in time, so set time dimension to 1
    field_block_[fieldID][3] = 1;


    // dimensions of 4D field data in file
    // each proc will write a field block, so need field_block*nProc in each dimension
    filespace_dims[0] = filespace_maxdims[0] = field_block_[fieldID][0]*nProcxyz_[0];
    filespace_dims[1] = filespace_maxdims[1] = field_block_[fieldID][1]*nProcxyz_[1];
    filespace_dims[2] = filespace_maxdims[2] = field_block_[fieldID][2]*nProcxyz_[2];
    // initialize time dimension to totWrites_+1
    filespace_dims[3] = totWrites_ ;
    // let time dimension be extensible (for restarts?)
    filespace_maxdims[3] = H5S_UNLIMITED;

    // create file dataspace for this field
    filespace_[fieldID] = H5Screate_simple(ndims_, filespace_dims, filespace_maxdims);
    assert(filespace_[fieldID]>=0);

    // set up file dataspace to be chunked
    hid_t create_chunks_plist = H5Pcreate(H5P_DATASET_CREATE);
    // make chunks the size of the field blocks
    H5Pset_chunk(create_chunks_plist, ndims_, field_block_[fieldID]);
    assert(create_chunks_plist>=0);

    // determine label for each field
    switch(fieldID) {
	case 0: strcpy(fieldname,"Ex"); break;
	case 1: strcpy(fieldname,"Ey"); break;
	case 2: strcpy(fieldname,"Ez"); break;
	case 3: strcpy(fieldname,"Bx"); break;
	case 4: strcpy(fieldname,"By"); break;
	case 5: strcpy(fieldname,"Bz"); break;
	case 6: strcpy(fieldname,"Jx"); break;
	case 7: strcpy(fieldname,"Jy"); break;
	case 8: strcpy(fieldname,"Jz"); break;
	case 9: strcpy(fieldname,"Bx_tm1"); break;
	case 10: strcpy(fieldname,"By_tm1"); break; 
	case 11: strcpy(fieldname,"Bz_tm1"); break;
	case 12: strcpy(fieldname,"rho"); break;
    }
    // create a dataset for each field
    field_dataset_[fieldID] = H5Dcreate(fields_group_id_, fieldname, H5T_NATIVE_DOUBLE, 
                               filespace_[fieldID], H5P_DEFAULT, create_chunks_plist, H5P_DEFAULT);
    assert(field_dataset_[fieldID]>=0);
    H5Pclose(create_chunks_plist); 

    // above we have prepared the hdf5 file for writing.
    // now we need to define the dataspace of each field in the program memory
    grid->getnxyzTot(nxyzTot);
    memspace_dims[0] = nxyzTot[0];
    memspace_dims[1] = nxyzTot[1];
    memspace_dims[2] = nxyzTot[2];
    memspace_dims[3] = 1;
    memspace_[fieldID] = H5Screate_simple(ndims_, memspace_dims, NULL);
    assert(memspace_[fieldID]>=0);
  }
  
  H5Gclose(fields_group_id_);
  delete [] phys_dims;
  delete [] filespace_dims;
  delete [] filespace_maxdims;
  delete [] nxyzTot;
  delete [] memspace_dims;
}

FieldTimeseriesIO::~FieldTimeseriesIO() {
  for(int i=0; i<nFieldDatasets_; i++) {
    H5Sclose(memspace_[i]);
    H5Sclose(filespace_[i]);
    H5Dclose(field_dataset_[i]);
    delete[] field_block_[i];
  }
  free(memspace_);
  free(filespace_);
  free(field_dataset_);
  delete[] field_block_;
}

int FieldTimeseriesIO::writeAField(const int fieldID, double*** field_data, const int iwrite) {
  assert(iwrite<=totWrites_);
  int status = -1;

  // Note: field_block_, filespace_, memspace_ allocated and set in constructor

  hsize_t* offset = new hsize_t[ndims_];
  hsize_t* stride = new hsize_t[ndims_];
  hsize_t* count = new hsize_t[ndims_];
  assert(offset!=NULL);
  assert(stride!=NULL);
  assert(count!=NULL);

  // select the subset of the file dataspace that this proc will be writing to
  offset[0] = myijk_[0]*field_block_[fieldID][0];
  offset[1] = myijk_[1]*field_block_[fieldID][1];
  offset[2] = myijk_[2]*field_block_[fieldID][2];
  offset[3] = iwrite;
  stride[0] = stride[1] = stride[2] = stride[3] = 1;
  count[0] = count[1] = count[2] = count[3] = 1;
  
  status = H5Sselect_hyperslab(filespace_[fieldID], H5S_SELECT_SET, 
				   offset, stride, count, field_block_[fieldID]);
  assert(status>=0);

  // select the subset of the memory dataspace that we are writing from
  offset[0] = offset[1] = offset[2] = 1;  // skip first cell in each spatial dim, which is ghost
  offset[3] = 0;
  stride[0] = stride[1] = stride[2] = stride[3] = 1;
  count[0] = count[1] = count[2] = count[3] = 1;
  
  status = H5Sselect_hyperslab(memspace_[fieldID], H5S_SELECT_SET, 
				   offset, stride, count, field_block_[fieldID]);
  assert(status>=0);

  // write from field_data in memspace to filespace in file
  status = H5Dwrite(field_dataset_[fieldID], H5T_NATIVE_DOUBLE, 
		      memspace_[fieldID], filespace_[fieldID],
		      data_xfer_plist_, field_data); 
  assert(status>=0);

  delete [] offset;
  delete [] stride;
  delete [] count;
  
  return status;
}

int FieldTimeseriesIO::writeFields(Grid* grid, const int which, const int iwrite){

  double**** fieldPtr = grid->getFieldPtr();  
  assert(fieldPtr!=NULL);

  //if(input->write_all_fields || input->write_E) {
  if(which==4 || which==1) {
    writeAField(grid->getExID(), fieldPtr[grid->getExID()], iwrite);
    writeAField(grid->getEyID(), fieldPtr[grid->getEyID()], iwrite);
    writeAField(grid->getEzID(), fieldPtr[grid->getEzID()], iwrite);
  }
  //if(input->write_all_fields || input->write_B) {
  if(which==4 || which==2) {
    writeAField(grid->getBxID(), fieldPtr[grid->getBxID()], iwrite);
    writeAField(grid->getByID(), fieldPtr[grid->getByID()], iwrite);
    writeAField(grid->getBzID(), fieldPtr[grid->getBzID()], iwrite);
  }
  //if(input->write_all_fields || input->write_J) {
  if(which==4 || which==3) {
    writeAField(grid->getJxID(), fieldPtr[grid->getJxID()], iwrite);
    writeAField(grid->getJyID(), fieldPtr[grid->getJyID()], iwrite);
    writeAField(grid->getJzID(), fieldPtr[grid->getJzID()], iwrite);
  }
  //if(input->write_all_fields || input->write_rho) {
  if(which==4 || which==0) {
    writeAField(grid->getrhoID(), fieldPtr[grid->getrhoID()], iwrite);
  }

  return 0;
}



