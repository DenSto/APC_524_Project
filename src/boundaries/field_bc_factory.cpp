#include "../globals.hpp"
#include "field_bc_factory.hpp"

/*! Construct the boundary condition array (must be freed!)
 * Takes in an array of size 6. */
void Field_BC_Factory::Construct(Domain* domain, Grid *grids, Input_Info_t *input_info){

    const char (*bound)[NCHAR] = input_info->fields_bound; 
    // convert c string for MPI to std::string for BC_Particle
    // size NCHAR is in correspondence definition in Input_Info_t
    std::string types[6];
    for(int i=0;i<6;i++){
        types[i].assign(bound[i]);
        //std::cerr << types[i] <<std::endl;
    }

    BC_Field** ret = new BC_Field*[6];
    assert(ret != NULL);
#if USE_MPI
    std::string periodic ("periodic"); //periodic gets treated specially if using MPI
    std::string mpi ("MPI");
    int* nProc = domain->getnProcxyz();	
    int* myLoc = domain->getmyijk();	
#endif

    for(int i = 0; i < 3; i++){
#if USE_MPI
        if(nProc[i] > 1){ // needs MPI (i.e. periodic BCs will be treated by MPI)
	    int partitionIndex = myLoc[i];
	    short inMiddle = (partitionIndex != 0) && (partitionIndex != nProc[i] - 1);
            short onRight = (partitionIndex == nProc[i]-1); // used to break MPI deadlock
            if(debug)fprintf(stderr,"rank=%d:dim=%d,onRight=%d\n",rank_MPI,i,onRight);
	    // compare return 0 when equal
	    if(periodic.compare(types[2*i])==0 || inMiddle){ // Two MPI communication loops
		// left boundary condition first, unless onRight
		ret[2*i+onRight]=lookup(mpi)(-(i+1),domain,grids,input_info);
		// Right boundary condition first, unless onRight
		ret[2*i+1-onRight]=lookup(mpi)(i+1,domain,grids,input_info);
            } else { // One MPI communication loop
		int MPIisRight;	
		// Physical side (to be calculated first to avoid MPI deadlock !!)
		if(partitionIndex == 0){ // Physical boundary on left
                    ret[2*i]=lookup(types[2*i])(-(i+1),domain,grids,input_info);
		    MPIisRight=1;
		} else { // Physical boundary on right
		    ret[2*i]=lookup(types[2*i+1])(i+1,domain,grids,input_info);
		    MPIisRight=-1;
		}
		// MPI side (left or right, compute second for efficiency)
		ret[2*i + 1]=lookup(mpi)(MPIisRight*(i+1),domain,grids,input_info);
	    }
 	} else // treat serially
#endif
	{
	    // Left boundary condition
	    ret[2*i]=lookup(types[2*i])(-(i+1),domain,grids,input_info);
            // Right boundary condition
	    ret[2*i+1]=lookup(types[2*i+1])(i+1,domain,grids,input_info);
	}
    }			 
    
    grids->setBoundaries(ret);    

}
