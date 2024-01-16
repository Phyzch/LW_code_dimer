#include"util.h"
#include"system.h"
using namespace std;

#define pi2 3.141592653589793*2 

//// allocate space for matrix and array
system::system() {
;
}



system::~system(){
    // destructor for system()
    delete [] x_electronic;
    delete [] y_electronic;
    delete [] electronic_state_energy;
}