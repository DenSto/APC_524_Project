#include "particle_bc_factory.hpp"
   /*
   	* Construct the boundary condition array (must be freed!)
	* Takes in an array of size 6.
    */
     //BC_Particle** constructConditions(Domain* domain, const std::string* types){
     BC_Particle** Part_BC_Factory::constructConditions(Domain* domain, Input_Info_t* info){
                // convert c string for MPI to std::string for BC_Particle
                // size NCHAR is in correspondence definition in Input_Info_t

		setInfo(info);
		std::string types[6];
		for(int i=0;i<6;i++){
		    types[i].assign(info->parts_bound[i]);
		    //std::cerr << types[i] <<std::endl;
		}
#if USE_MPI
		std::string periodic ("periodic"); //periodic gets treated specially if using MPI
		std::string mpi ("MPI");
		int* nProc = domain->getnProcxyz();	
		int* myLoc = domain->getmyijk();	
#endif
		BC_Particle** ret = new BC_Particle*[6];
		assert(ret != NULL);
		for(int i = 0; i < 3; i++){
#if USE_MPI
			if(nProc[i] > 1){ // needs MPI (i.e. periodic BCs will be treated by MPI)
				int partitionIndex = myLoc[i];
				short inMiddle = (partitionIndex != 0) && (partitionIndex != nProc[i] - 1);
				// compare return 0 when equal
				if(periodic.compare(types[2*i])==0 || inMiddle){ // Two MPI communication loops
					// both boundaries are handled in each ret, to avoid MPI deadlock
					// Left boundary condition
					ret[2*i]=lookup(mpi)(domain,i,0,types[2*i]);

					// Right boundary condition
					ret[2*i+1]=lookup(mpi)(domain,i,1,types[2*i+1]);
				} else { // One MPI communication loop
					// Physical side (to be calculated first!!)
					int MPIisRight = 0;	
					if(partitionIndex == 0){ // Physical boundary on left
						ret[2*i]=lookup(types[2*i])(domain,i,0,types[2*i]);
						MPIisRight=1;
					} else { // Physical boundary on right
						ret[2*i]=lookup(types[2*i+1])(domain,i,1,types[2*i+1]);
						MPIisRight=0;
					}
					// MPI side (left or right, compute second for efficiency)
					ret[2*i + 1]=lookup(mpi)(domain,i,MPIisRight, types[2*i + MPIisRight]);
				}
 			} else // treat serially
#endif
			{
				// Left boundary condition
				ret[2*i]=lookup(types[2*i])(domain,i,0,types[2*i]);

				// Right boundary condition
				ret[2*i+1]=lookup(types[2*i+1])(domain,i,1,types[2*i+1]);
			}
		}			 
		return ret;
	}
