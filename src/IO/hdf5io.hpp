#ifndef HDF5IO_HPP
#define HDF5IO_HPP

#if USE_MPI
#include "mpi.h"
#endif

#include "../domain/domain.hpp"
#include "../grid/grid.hpp"
#include "input.hpp"

#include "hdf5.h"

/*! Class that creates an hdf5 file and sets basic props */
class Hdf5IO {
    public:
	Hdf5IO(void) {};
        Hdf5IO(const char* filename);
        ~Hdf5IO(void);

	hid_t getFileID() 
	  { return file_id_; }
	hid_t getFileAccessPlist() 
	  {return file_access_plist_;}
	hid_t getDataXferPlist()
	  {return data_xfer_plist_;}

    protected:
	// file id
	hid_t file_id_;
	
	// property list ids
	hid_t file_access_plist_;
	hid_t data_xfer_plist_;

};

/*! Class for writing fields in time series to hdf5 */
class FieldTimeseriesIO : public Hdf5IO {
    public: 
	FieldTimeseriesIO(Hdf5IO* io, Grid* grid, Domain* domain, Input_Info_t* input, const int totWrites);
	~FieldTimeseriesIO(void);
	int getnFieldDatasets() {return nFieldDatasets_;}
	int writeAField(const int fieldID, double*** data, const int iwrite);
	int writeFields(Grid* grid, Input_Info_t* input, const int iwrite);

    private:
	// number of field datasets
	const int nFieldDatasets_;

	// total number of times fields will be written
	const int totWrites_;

	// number of dimensions
	const int ndims_;

	// information about proc layout
	// number of procs in each dimension
	int* nProcxyz_;

	// where on proc grid 
	int* myijk_;

	// group id
	hid_t fields_group_id_;

	// memory dataspace ids for fields
	hid_t* memspace_;

	// file dataspace ids for fields
	hid_t* filespace_;

        // dataset ids
	hid_t* field_dataset_;

	// fields written by blocks
	// each field has its own block dimensions,
	// given by field_block_[fieldID]
	hsize_t** field_block_;

};

#endif
