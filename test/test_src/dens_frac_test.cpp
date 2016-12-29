#include <iostream>
using namespace std;

int main(){
    int nspec;
    cout << "Eneter number of species\n";
    cin >> nspec;

    cout << "Enter fractional densities\n";
    double *dens = new double[nspec];
    double cden = 0.0; // cummulative density
    for(int i=0;i<nspec;i++){
        cin >> dens[i];
        cden += dens[i];
    }

    cout << "The normalized fractional densities are " << endl;
    for(int i=0;i<nspec;i++){
        dens[i]/=cden;
        cout << "dens[" << i << "]=" << dens[i] << endl;
    }

    int npart;
    cout << "Enter total number of particles\n";
    cin >> npart;

    int ispec = 0;
    cden = dens[0];
    for(int ip=0;ip<npart;ip++){
        if(ip >= cden*npart){
            ispec += 1;
            cden  += dens[ispec];
	}
        cout << ip << "th particle is of type " << ispec << endl;
    }

    return 0;
}
