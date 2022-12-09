//
// Created by phyzch on 4/14/20.
//Used for putting non essential function in same file.
//
#include"system.h"
#include"util.h"
using namespace std;
namespace fs = std::experimental::filesystem;

// check some dimension parameter
void full_system::dimension_check() {

    if (s.electronic_state_num == 1) {
        output << "Global Matrix: 2*" << d.dmatsize[0] << " = " << matsize << endl;
    }
    else if (s.electronic_state_num == 2) {
        output << "Global Matrix: " << total_matsize << endl;
    }
    output << "off-diagonal matrix number  " << total_offnum << endl;
    output << "Whole matrix element number  " << total_matnum << endl;
    log << "off-diagonal matrix number  " << total_offnum << endl;
    log << "Whole matrix element number  " << total_matnum << endl;
};


full_system::~full_system(){
    // release the space allocated by new. destructor.
    int i;
    delete [] matsize_each_process;
    delete [] mat_offnum_each_process;
    delete [] matnum_each_process;
    delete [] matsize_offset_each_process;
    delete [] matnum_offset_each_process;

    delete [] dstate;
    delete [] dstate_all;

    delete [] remoteVecCount;
    delete [] remoteVecPtr;
    delete [] remoteVecIndex;
    delete [] tosendVecCount;
    delete [] tosendVecPtr;
    delete [] tosendVecIndex;
    delete [] recv_x;
    delete [] recv_y;
    delete [] send_x;
    delete [] send_y;

}

template <typename T > void broadcast_1d_vector(vector<T> & vector_array, int & array_size, MPI_Datatype datatype , int root){
    int i;

    MPI_Bcast(&array_size, 1, MPI_INT, root, MPI_COMM_WORLD);

    T *array = new T [array_size];
    if (my_id == root){
        for(i=0;i<array_size;i++) {
            array[i] = vector_array[i];
        }
    }
    MPI_Bcast(&array[0], array_size, datatype, root, MPI_COMM_WORLD);
    if (my_id != root){
        for(i=0;i<array_size;i++){
            vector_array.push_back(array[i]);
        }
    }

    delete [] array;
}
// explicitly instantiate template.
// see :https://bytefreaks.net/programming-2/c/c-undefined-reference-to-templated-class-function
template void broadcast_1d_vector <int>( vector<int> & vector_array, int & array_size, MPI_Datatype datatype, int root);
template void broadcast_1d_vector <double> (vector<double> & vector_array , int & array_size, MPI_Datatype datatype, int root);



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