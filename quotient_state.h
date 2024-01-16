#include<iostream>
#include"util.h"//
// Created by phyzch on 6/22/20.
//

#ifndef QUANTUM_MEASUREMENT_QUOTIENT_STATE_H
#define QUANTUM_MEASUREMENT_QUOTIENT_STATE_H
using namespace std;


struct quotient_state {   // monomer quotient_space_state.

    // each quotient_state is defined for pair (exciton_state_index_list, vmode), here vmode is vibrational mode in another monomer
    // exciton_state_index_list is exciton state ( 0 or 1).
    // states (exciton_state_index_list, vmode1, vmode2) is grouped into different group according to (exciton_state_index_list, vmode2) for monomer1_quotient_state_list and (exciton_state_index_list, vmode1) for monomer2_quotient_state_list
    // monomer1_quotient_state_list is used when we construct anharmonic coupling in monomer1, states are grouped according to (exciton_state_index_list, vmode2)
    // monomer2_quotient_state_list is used when we construct anharmonic coupling in monomer2, states are grouped according to (exciton_state_index_list, vmode1)

    // Take monomer1_quotient_state_list for example:
    int exciton_state;  // exciton state. denote different potential energy surface
    vector<int> vmode;  // vibrational mode of monomer2


    vector <int> full_hamiltonian_state_index_list;  // list of state index in exciton_state_index_list + monomer vib state wave function that is defined with (exciton_state_index_list, vmode1, vmode2),
                                                // which is grouped according to (exciton_state_index_list, vmode2)

    vector <int> monomer_state_index_list; // sorted list. records index in monomer1 reduced density matrix basis. (monomer states with vibrational quantum number vmode1)


    // anharmonic_coupling_info_index records anharmonic coupling for vib states in monomer 1.
    // list of tuple (i,j,k,l,m):
    // i : vib state in monomer1, j : vib state in monomer1.  i,j monoer state coupled with each other anharmonically.
    // k: state index in full_matrix, l: state index in full matrix.
    // m : index in monomer Hamiltonian for local anharmonic coupling
    vector<vector<int>> anharmonic_coupling_info_index_list;

    vector<double> anharmonic_coupling_value_list;

    // initialize quotient state. defined with vibrational states in another monomer and exciton states
    quotient_state(vector<int>  & vmode1, int exciton_state1){
        exciton_state = exciton_state1;
        vmode = vmode1;
    }
};

int binary_insert_dxindex(vector <int> & list, int key);
int compare_quotient_state(int sys_state1, const vector <int>  & vmode1, int sys_state2, const vector <int> & vmode2);
int binary_search_monomer_state_index(const vector <int> & list, int key, bool & exist);
void insert_quotient_state(vector <quotient_state> & quotient_states_list, int exciton_state, vector<int> & vmode1, int full_hamiltonian_state_index, int monomer_state_index);
void save_detector_quotient_state_data_for2( const vector <quotient_state> & d1list, const vector <quotient_state> & d2list, string path);
void load_detector_quotient_state_data_for2( vector<quotient_state> & d1list, vector <quotient_state> & d2list, string path);




int find_position_for_insert_binary(const vector<vector<int>> & vmode, const vector<int> & ndetector, bool & exist);  // write in compute_matrix_energy_window.cpp
int find_location_binarysearch_quotient_state(const vector<quotient_state> & list, int sys_state, const vector <int> & vmode, bool & exist);

#endif //QUANTUM_MEASUREMENT_QUOTIENT_STATE_H
