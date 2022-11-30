// Chenghao Zhang  cz38@illinois.edu
#include"system.h"
#include"util.h"
using namespace std;
namespace fs = std::experimental::filesystem;

// Information for MPI program
int my_id;
int num_proc;

int main(int argc,char * argv []) {
    srand(time(0));
    string parentpath= "/home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/BChl_dimer_model/BChl try";

    string path;

    string s;
    string s1;
    int i;
    int Filenumber=1;

    // MPI Command
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,&num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);

    for(i=0;i<Filenumber;i++){
        path=parentpath;

        { // the parenthese here let destructor called after we use this instance.
            full_system photon_entangled_system(path);  // set parameter and construct Hamiltonian.
            photon_entangled_system.Quantum_evolution(); // creat initial state (or read from file). Then complete simulation.
        }

    }

    MPI_Finalize();
}


