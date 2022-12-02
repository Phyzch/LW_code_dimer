#include"system.h"
#include"util.h"
using namespace std;
// add comment to test git.

// if we are going to use energy_window, we don't have to pay attention to detdim and dmatdim, because we will calculate these number in our own functioin
// however, if we are not going to add energy window, we have to make sure these values are right before begin our simulation.

detector::detector(){
    ; // do nothing here, I write this code because I don't know what will compiler do if I don't initialize class.
}

detector::~detector(){
    int i,m,j,k;
    int max_qn;

    delete [] deln;
    delete [] nbar;
    delete [] dirow;
    delete [] dicol;
    delete [] dv;
    delete [] xd;
    delete [] yd;
    delete [] proptime;
    delete [] initial_state_energy;

    for(i=0; i < electronic_state_num; i++){
        delete [] dmatsize_each_process[i];
        delete [] doffnum_each_process[i];
        delete [] dmatnum_each_process[i];
        delete [] dmat_offset_each_process[i];
        delete [] total_dmat[i];
        delete [] total_dirow[i];
        delete [] total_dicol[i];
        delete [] xd_all[i];
        delete [] yd_all[i];

        delete [] mfreq[i];
        delete [] electron_phonon_coupling[i];
        delete [] initial_detector_state[i];
        delete [] aij[i];

    }
    delete [] initial_state_pc_id;
    delete [] initial_state_index;
    delete [] initial_detector_state;

    delete [] aij;
    delete [] mfreq;
    delete [] electron_phonon_coupling;
    delete [] dmatsize_each_process;
    delete [] doffnum_each_process;
    delete [] dmatnum_each_process;
    delete [] dmat_offset_each_process;
    delete [] total_dmat;
    delete [] total_dirow;
    delete [] total_dicol;

    delete [] xd_all;
    delete [] yd_all;

    // ------ delete franck condon factor table
    for(i=0;i<nmodes[0];i++){
        if (max_qn < nmax[0][i]) {
            max_qn = nmax[0][i];
        }
    }
    max_qn = max_qn + 1;
    for(m=0; m < electronic_state_num;m++){
        for(j=0;j<nmodes[m];j++){
            for(k=0;k<max_qn;k++){
                delete [] franck_condon_factor_table[m][j][k];
            }
            delete [] franck_condon_factor_table[m][j];
        }
        delete [] franck_condon_factor_table[m];
    }
    delete [] franck_condon_factor_table;
    // -----------------

    for(m=0;m<electronic_state_num;m++){
        delete[] nmax[i];
    }
    delete [] nmax;
    delete [] dmatsize;
    delete [] dmatnum;
    delete [] doffnum;
    delete [] nmodes;
}
