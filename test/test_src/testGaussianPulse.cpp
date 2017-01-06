#include <iostream>
#include "GaussianPulse.hpp"

using namespace std;

int main(){
    int nwaves;
    cout << "Eneter number of waves\n";
    cin >> nwaves;

    double *amp = new double[nwaves];
    cout << "Enter amplitudes\n";
    for(int i=0;i<nwaves;i++){
        cin >> amp[i];
    }

    double *omg= new double[nwaves];
    cout << "Enter omegas\n";
    for(int i=0;i<nwaves;i++){
        cin >> omg[i];
    }

    double *pha = new double[nwaves];
    cout << "Enter phases in degree\n";
    for(int i=0;i<nwaves;i++){
        cin >> pha[i];
    }

    double *del = new double[nwaves];
    cout << "Enter delyas\n";
    for(int i=0;i<nwaves;i++){
        cin >> del[i];
    }

    double *dwt = new double[nwaves];
    cout << "Enter widths\n";
    for(int i=0;i<nwaves;i++){
        cin >> dwt[i];
    }

    Gaussian_Pulses_t *pulses = new Gaussian_Pulses_t;
    pulses->nwaves   = nwaves;
    pulses->peakamps = amp;
    pulses->omegas   = omg;
    pulses->phases   = pha;
    pulses->delays   = del;
    pulses->invWidths= dwt; 

    double t=0.0;
    do{
        cout << "Enter time to compute. To abort, enter t<-999.9\n";
        cin >> t;
        cout << "pulse value = " << GaussianPulses(t,pulses) << endl; 
    }while(t>-999.9);

    return 0;
}
