// Chenghao Zhang  cz38@illinois.edu
#include"system.h"
#include"util.h"
using namespace std;
namespace fs = std::experimental::filesystem;
void check_and_create_file(string parent_path, string path);

// Information for MPI program
int my_id;
int num_proc;

// About matflag in input.txt: If matflag==2, +We output all x,y (after the pre_coupling), matrix element, detector matrix element etc.
// if matflag==1: We don't output anything but still we will save our final simulation results in save.txt
int main(int argc,char * argv []) {
    srand(time(0));
    string parentpath= "/home/phyzch/CLionProjects/Quantum_Measurement_MPI_version/result/CVPT project/"
                       "Tunneling between molecules/two side different dof/Batch simulation random frequency/state 0/c=0.1/V=80/";
    string cvpt_parent_path = "/home/phyzch/CLionProjects/CVPT/data/Potential for flow between molecules/"
                              "batch simulation no self-coupling from potential/0.2/";
//    string cvpt_parent_path = "/home/phyzch/CLionProjects/CVPT/data/New Criteria/cutoff 0.05 GOE/";
    string cvpt_path;
    string path;
    string path1;
    string path2;
    string path3;
    string mode_path ;
    string s;
    string s1;
    int i,j,k,l;
    int Filenumber=1;

    // MPI Command
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,&num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);

    for(i=0;i<Filenumber;i++){
        path=parentpath;
        cvpt_path =cvpt_parent_path;

        { // the parenthese here let destructor called after we use this instance.
            full_system photon_entangled_system(path,cvpt_path);  // set parameter and construct Hamiltonian.
            photon_entangled_system.Quantum_evolution(); // creat initial state (or read from file). Then complete simulation.
        }

    }

    MPI_Finalize();
}

void get_current_path(){
    char buff[FILENAME_MAX];
    getcwd(buff,FILENAME_MAX);
    string current_location(buff);
    cout<<"Current location is :"<<current_location<<endl;   // code for debug output the current location
}

void check_and_create_file(string parent_path, string path){
    // check existence of sub-folder, if not , create it.
    // copy input.txt into subfolder to do simulation.
    struct stat statbuf;
    bool isDir = false;
    if (stat(path.c_str(), &statbuf) != -1) {
        // get permission to access directory
        if (S_ISDIR(statbuf.st_mode)) {
            // is directory
            isDir = true;
        }
    }
    if (isDir){
        ;  // directory already exists
    }
    else {// create directory
        if (!mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)) {
            printf("File created!");  // successfully create folder.
        } else {
            printf("Fail to create the file, directory may not exist.");
            exit(-1);
        }
    }
    // create/update input.txt from existing folder to subfolder for simulation.
    const fs::path dst= path + "input.txt";
    const fs::path src= parent_path + "input.txt";
    if(fs::exists(src)) {
        // copy input.txt into subfolder.
        experimental::filesystem::copy_file(src, dst,fs::copy_options::overwrite_existing);
    }
    else{
        cout<<"File do not exist."<<endl;
        exit(-1);
    }
}

