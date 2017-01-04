#include "field_bc_factory.hpp"

/*! Construct the boundary condition array (must be freed!)
 * Takes in an array of size 6. */
void Field_BC_Factory::constructConditions(Domain* domain, Grid *grids, const char (*bound)[NCHAR]){
    // convert c string for MPI to std::string for BC_Particle
    // size NCHAR is in correspondence definition in Input_Info_t
    std::string types[6];
    for(int i=0;i<6;i++){
        types[i].assign(bound[i]);
        //std::cerr << types[i] <<std::endl;
    }
#if USE_MPI
    std::string periodic ("periodic"); //periodic gets treated specially if using MPI
    std::string mpi ("MPI");
    int* nProc = domain->getnProcxyz();	
    int* myLoc = domain->getmyijk();	
#endif
    BC_Field** ret = new BC_Field*[6];
    assert(ret != NULL);
    for(int i = 0; i < 3; i++){
#if USE_MPI
        if(nProc[i] > 1){ // needs MPI (i.e. periodic BCs will be treated by MPI)
	    int partitionIndex = myLoc[i];
	    short inMiddle = (partitionIndex != 0) && (partitionIndex != nProc[i] - 1);
	    // compare return 0 when equal
	    if(periodic.compare(types[2*i])==0 || inMiddle){ // Two MPI communication loops
		// Left boundary condition
		ret[2*i]=lookup(mpi)(domain,grids,-(i+1));
		// Right boundary condition
		ret[2*i+1]=lookup(mpi)(domain,grids,i+1);
            } else { // One MPI communication loop
		int MPIisRight;	
		// Physical side (to be calculated first!!)
		if(partitionIndex == 0){ // Physical boundary on left
                    ret[2*i]=lookup(types[2*i])(domain,grids,-(i+1));
		    MPIisRight=1;
		} else { // Physical boundary on right
		    ret[2*i]=lookup(types[2*i+1])(domain,grids,i+1);
		    MPIisRight=-1;
		}
		// MPI side (left or right, compute second for efficiency)
		ret[2*i + 1]=lookup(mpi)(domain,grids,MPIisRight*(i+1));
	    }
 	} else // treat serially
#endif
	{
	    // Left boundary condition
	    ret[2*i]=lookup(types[2*i])(domain,grids,-(i+1));
            // Right boundary condition
	    ret[2*i+1]=lookup(types[2*i+1])(domain,grids,(i+1));
	}
    }			 
    
    grids->setBoundaries(ret);    

}
