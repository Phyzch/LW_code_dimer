#include"system.h"
#include"util.h"
using namespace std;
// add comment to test git.

// if we are going to use energy_window, we don't have to pay attention to detdim and dmatdim, because we will calculate these number in our own functioin
// however, if we are not going to add energy window, we have to make sure these values are right before begin our simulation.

monomer::monomer(){
    ; // do nothing here, I write this code because I don't know what will compiler do if I don't initialize class.
}

monomer::~monomer(){
    int i,m,j,k;
    int max_qn;

    delete [] deln;
    delete [] nbar;
    delete [] monomer_irow;
    delete [] monomer_icol;
    delete [] monomer_vibrational_states_quantum_number_list;
    delete [] xd;
    delete [] yd;
    delete [] proptime;
    delete [] initial_state_energy;

    for(i=0; i < electronic_state_num; i++){
        delete [] monomer_matsize_each_process[i];
        delete [] monomer_matsize_offset_each_process[i];
        delete [] monomer_offnum_each_process[i];
        delete [] monomer_matnum_each_process[i];
        delete [] monomer_mat_offset_each_process[i];
        delete [] total_monomer_mat[i];
        delete [] total_monomer_irow[i];
        delete [] total_monomer_icol[i];
        delete [] xd_all[i];
        delete [] yd_all[i];

        delete [] mfreq[i];
        delete [] electron_phonon_coupling[i];
        delete [] initial_vibrational_state[i];
        delete [] aij[i];

    }
    delete [] initial_state_pc_id;
    delete [] initial_state_index;
    delete [] initial_vibrational_state;

    delete [] aij;
    delete [] mfreq;
    delete [] electron_phonon_coupling;
    delete [] monomer_matsize_each_process;
    delete [] monomer_matsize_offset_each_process;
    delete [] monomer_offnum_each_process;
    delete [] monomer_matnum_each_process;
    delete [] monomer_mat_offset_each_process;
    delete [] total_monomer_mat;
    delete [] total_monomer_irow;
    delete [] total_monomer_icol;

    delete [] xd_all;
    delete [] yd_all;

    for(m=0;m<electronic_state_num;m++){
        delete[] nmax[m];
    }
    delete [] nmax;
    delete [] monomer_matsize;
    delete [] monomer_matnum;
    delete [] monomer_offnum;
    delete [] nmodes;
}
