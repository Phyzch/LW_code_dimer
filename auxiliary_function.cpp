//
// Created by phyzch on 4/14/20.
//Used for putting non essential function in same file.
//
#include"system.h"
#include"util.h"
using namespace std;

// check some dimension parameter
void full_system::dimension_check() {
    double errflag = 0;
    if (errflag != 0) {
        log << " Dimension Problem: Error flag=" << errflag << endl;
        cout<<"Dimension error, Error flag="<<errflag;
        exit(-1);
    }

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