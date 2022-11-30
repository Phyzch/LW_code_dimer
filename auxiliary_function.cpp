//
// Created by phyzch on 4/14/20.
//Used for putting non essential function in same file.
//
#include"system.h"
#include"util.h"
using namespace std;

// code used to reform the final output if we re-start our simulation
void full_system::replace_first_line() {
    // This code is quite dumb... I find if I continue my simulation, I have to correct the endtime I give. so I have to rewrite whole output.txt for the first line.
    ifstream fin;
    ofstream temp;
    int x;
    if(my_id==0) {
        string old_path = path + "output.txt";
        fin.open(old_path);
        string new_path = path + "outputnew.txt";
        temp.open(new_path);
        string line;
        getline(fin, line); // read first line
        temp << delt << " "  << tmax << " " << tprint << endl;
        while (std::getline(fin, line)) {
            if (line != "") {
                temp << line << endl;
            }
        }
        temp.close();
        fin.close();
        if (remove(old_path.c_str()) == 0) {
            rename(new_path.c_str(), old_path.c_str());
        } else {
            cout << "Remove file failed" << endl;
            cin >> x;
        }
    }
}

// check some dimension parameter
void full_system::dimension_check() {
    double errflag = 0;
    if (errflag != 0) {
        log << " Dimension Problem: Error flag=" << errflag << endl;
        cout<<"Dimension error, Error flag="<<errflag;
        exit(-1);
    }

    if (s.tlnum == 1) {
        output << "Global Matrix: 2*" << d.dmatsize[0] << " = " << matsize << endl;
    }
    else if (s.tlnum == 2) {
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
    delete [] sdnum;
    delete [] dstate;
    delete [] sdindex;
    delete [] sdmode;
    delete [] remoteVecCount;
    delete [] remoteVecPtr;
    delete [] remoteVecIndex;
    delete [] tosendVecCount;
    delete [] tosendVecPtr;
    delete [] tosendVecIndex;

    delete [] remote_vec_count_for_detenergy;
    delete [] remote_vec_ptr_for_detenergy;
    delete [] remote_vec_index_for_detenergy;
    delete [] to_send_vec_count_for_detenergy;
    delete [] to_send_vec_ptr_for_detenergy;
    delete [] to_send_vec_index_for_detenergy;
    delete [] x_for_detenergy;
    delete [] y_for_detenergy;
    delete [] send_x_for_detenergy;
    delete [] send_y_for_detenergy;
}