#ifndef HDF5IO_HPP
#define HDF5IO_HPP

#if USE_MPI
#include "mpi.h"
#endif

#include "../domain/domain.hpp"
#include "../grid/grid.hpp"
#include "input.hpp"

#include "hdf5.h"

enum fields_to_write {
  rho = 0,
  E = 1,
  J = 2,
  B = 3,
  ALL = 4
};

class FieldTimeseriesIO;

/*! Class that creates an hdf5 file and sets basic props */
class Hdf5IO {
    public:
        Hdf5IO(const char* filename, Grid* grid, Domain* domain, const int which_fields);
        ~Hdf5IO(void);

	hid_t getFileID() 
	  { return file_id_; }
	hid_t getFileAccessPlist() 
	  {return file_access_plist_;}
	hid_t getDataXferPlist()
	  {return data_xfer_plist_;}
	hid_t getFieldsGroupID() 
	  { return fields_group_id_; }
	int* getnProcxyz() 
	  { return nProcxyz_;}
	int* getmyijk()
	  { return myijk_;}

	int writeFields(Grid* grid, double time);
	int writeTime(double time);

    protected:
	// file id
	hid_t file_id_;
	
	// property list ids
	hid_t file_access_plist_;
	hid_t data_xfer_plist_;

	// time attribute
	hid_t t_dataset_;
	hid_t t_dataspace_;
	hid_t t_memspace_;

	// for field diagnostics
	hid_t fields_group_id_;
	const int which_fields_;
	// information about proc layout
	// number of procs in each dimension
	int* nProcxyz_;

	// where on proc grid 
	int* myijk_;
	
	FieldTimeseriesIO* Ex_tsio_;
	FieldTimeseriesIO* Ey_tsio_;
	FieldTimeseriesIO* Ez_tsio_;
	FieldTimeseriesIO* Bx_tsio_;
	FieldTimeseriesIO* By_tsio_;
	FieldTimeseriesIO* Bz_tsio_;
	FieldTimeseriesIO* Jx_tsio_;
	FieldTimeseriesIO* Jy_tsio_;
	FieldTimeseriesIO* Jz_tsio_;
	FieldTimeseriesIO* rho_tsio_;

};


/*! Class for writing fields in time series to hdf5 */
class FieldTimeseriesIO  {
    public: 
	FieldTimeseriesIO(Hdf5IO* io, Grid* grid, const int fieldID, std::string fieldname);
	~FieldTimeseriesIO(void);
	int writeField(double*** data);

    private:
	// number of dimensions
	const int ndims_;

	// information about proc layout
	// number of procs in each dimension
	int* nProcxyz_;

	// where on proc grid 
	int* myijk_;

	// file id
	hid_t file_id_;
	
	// property list ids
	hid_t data_xfer_plist_;

	// group id
	hid_t fields_group_id_;

	// memory dataspace ids for fields
	hid_t memspace_;

	// file dataspace ids for fields
	hid_t filespace_;

        // dataset ids
	hid_t field_dataset_;

	// fields written by blocks
	// each field has its own block dimensions,
	// given by field_block_[fieldID]
	hsize_t* field_block_;

};

#endif
