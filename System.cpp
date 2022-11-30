#include"util.h"
#include"system.h"
using namespace std;

#define pi2 3.141592653589793*2 

//// allocate space for matrix and array
system::system() {
    // x_electronic: real part of system wave function
    // y_electronic: image part of system wave function
    // electronic_state_energy: energy level of system's eigen state, also diagonal part of our Hamiltonian matrix
    // tlrho: density matrix of system
    // tlirow: row index for matrix element in tlmat
    // tlicol: column index for matrix element in tlmat
	x_electronic = new double[int(pow(2, tldim))];
    y_electronic = new double[int(pow(2, tldim))];
    electronic_state_energy = new double[tldim];
	tlmat = new double[int(pow(2, tldim))];
}



system::~system(){
    // destructor for system()
    delete [] x_electronic;
    delete [] y_electronic;
    delete [] electronic_state_energy;
    delete [] tlmat;
}