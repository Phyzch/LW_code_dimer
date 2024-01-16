//
// Created by phyzch on 6/19/20.
// Implement code to speed up calculation of reduced density matrix.
// Ideally if we only compute energy of monomer, the cost will be O( matsize)
//
#include "system.h"
#include"util.h"
#include <iostream>
using namespace std;


// we will construct vector of struct. use binary search to serch quotient_state q, if can't find one, insert it in vector list and insert index, dindex.
// if we find it already exist, still insert index , dindex.

// when store dindex and index, we also have to make sure it is sorted. So we can use binary search to search it. result store in q_index.
// first time we compute reduced density matrix, we have to record all reduced density matrix element we need to compute, and tell
// every quotient_state instance (i,j) where i and j is dindex they have to compute.

int compare_quotient_state(int sys_state1, const vector <int>  & vmode1, int sys_state2, const vector <int> & vmode2){
    // we compare vmode first. this way we don't have to merge sorting between different process.
    int size= vmode1.size();
    if(size != vmode2.size()){
        cout<<" Error, size of mode does not equal to each other!"<<endl;
        exit(-5);
    }
    int i;
    // now compare vmode.
    for(i=0;i<size;i++){
        if(vmode1[i]>vmode2[i]) return 1;
        else if(vmode1[i] < vmode2[i] ) return -1;
        else ;
    }
    if(sys_state1 > sys_state2) return 1;
    else if (sys_state1< sys_state2) return -1;
    return 0;
}

// method to use find_location_binary : create a list vector <quotient_state>.
// for example, to compute reduced density matrix for monomer 1, we construct quotient Hilbert space of Hilbert space H_{d1}
// record  x[state] and dx[state] corresponding to this quotient space state d1. (in full_hamiltonian_state_index_list and monomer_state_index_list).
int find_location_binarysearch_quotient_state(const vector<quotient_state> & list, int sys_state, const vector <int> & vmode, bool & exist){
    // used to sort quotient state list. return location to insert or manipulate, also return exist to indicate this quotient_state exist in list or not
    // use information of exciton_state_index_list and vmode to locate the state.
    // use binary search method.
    int right_flag= list.size()-1; // this step is crucial for binary search.
    int left_flag= 0;
    if(list.size()==0){
        exist= false;
        return 0;
    }
    int position= (left_flag + right_flag)/2;
    exist =false;
    int mark;
    while(left_flag < right_flag){
        mark = compare_quotient_state(sys_state, vmode, list[position].exciton_state, list[position].vmode);
        if(mark== 1){
            // state > list[position]
            left_flag=position+1;
            position= (left_flag + right_flag)/2;
        }
        else if (mark == -1){
            right_flag= position -1;
            position = (left_flag + right_flag )/2;
        }
        else if(mark==0){
            exist= true;
            return position;
        }
    }
    // Now situation is left_flag == right_flag == position or position == right_flag == left_flag -1. either case we compare list[position] with our state
    mark =  compare_quotient_state(sys_state, vmode, list[position].exciton_state, list[position].vmode);
    if( mark == 0 ){
        exist =true;
        return position; // return the position
    }
    else if (mark == 1){
        exist =false;
        return position + 1; // insert after the position
    }
    else{  // mark == -1.
        exist = false;
        return position;  // insert before the position.
    }
}

int compare_int(int a, int b){
    if (a>b) return 1;
    else if (a<b) return -1;
    else return 0;
}

int binary_search_monomer_state_index(const vector <int> & list, int key, bool & exist){
    // search key in list. If exist, set exist to true, else set to false.
    // return position of key if found.
    int left_flag=0;
    int right_flag = list.size()-1;
    exist=false;
    if(list.size()==0){
        cout<<"error! The quotient_state when initialized should have non-null full_hamiltonian_state_index_list and monomer_state_index_list" <<endl;
        exit(-7);
    }
    int position = (left_flag + right_flag)/2;
    int mark;
    while(left_flag<=right_flag){
        mark= compare_int(key,list[position]);
        if(mark == 1) {
            left_flag= position+1;
            position = (right_flag + left_flag)/2;
        }
        else if (mark ==-1){
            right_flag = position -1;
            position = (right_flag + left_flag)/2;
        }
        else{
            exist = true;
            return position;
        }
    }
    exist= false;
    return 0;
}

int binary_insert_dxindex(vector <int> & list, int key){
    // insert key into list.
    // we do not allow key already in list, if we find that, we report error.  (That's because for one quotient state we can't store same monomer state index and state index twice if program is right.)
    int left_flag=0;
    int right_flag = list.size()-1;
    if(list.size()==0){
        list.push_back(key);
        return 0;
    }
    int position = (left_flag + right_flag)/2;
    int mark;
    while(left_flag< right_flag){
        mark = compare_int(key,list[position]);
        if(mark == 1) {
            left_flag= position+1;
            position = (right_flag + left_flag)/2;
        }
        else if (mark ==-1){
            right_flag = position -1;
            position = (right_flag + left_flag)/2;
        }
        else{
            cout<<"Error, the index can not be the same when we construct quotient state."<<endl;
            exit(-6);
        }
    }
    mark = compare_int(key,list[position]);
    if (mark ==-1){
        list.insert(list.begin() + position, key);
        return position;
    }
    else if (mark==1){
        list.insert(list.begin() + position +1 ,key);
        return position +1;
    }
    else{
        cout<<"Error, the index can not be the same when we construct quotient state."<<endl;
        exit(-6);
    }
}

void insert_quotient_state(vector <quotient_state> & quotient_states_list, int exciton_state, vector<int> & vmode1, int full_hamiltonian_state_index, int monomer_state_index) {
    // monomer_state_index_list: index in monomer wavefunction for states in Hilbert space d.
    // full_hamiltonian_state_index_list: index in the whole system's wave function. (exciton state + monomer 1 vib state + monomer 2 vib state)

    // take quotient states for monomer1 for example.
    // vmode1: vibrational states in monomer2.  exciton_state_index_list: 0 / 1, representing different exciton state (different potential energy surface)


    bool exist;
    int quotient_state_position; // quotient_state_position to insert.
    int monomer_state_index_position; // quotient_state_position to insert full_hamiltonian_state_index_list, it should be at the same quotient_state_position

    // copy vmode1 to vmode
    vector<int> vmode;
    vmode = vmode1;

    // find whether the <quotient state> defined by pair (exciton_state_index_list, vmode) is already in the quotient_state_list
    quotient_state_position = find_location_binarysearch_quotient_state(quotient_states_list, exciton_state, vmode, exist); // find state marked by (exciton_state_index_list, vmode) if it exist in quotient_states_list.

    if(!exist){
        // quotient state does not exist before. We have to construct it and insert it into quotient_state quotient_states_list.
        quotient_state s (vmode, exciton_state); // construct quotient_state
        monomer_state_index_position = binary_insert_dxindex(s.monomer_state_index_list, monomer_state_index ); // insert monomer state index into quotient_states_list with order.
        s.full_hamiltonian_state_index_list.insert(s.full_hamiltonian_state_index_list.begin() + monomer_state_index_position, full_hamiltonian_state_index);
        // insert full hamiltonian state index into quotient_states_list at the same position as the monomer state index in  monomer_state_index_list.

        // insert quotient state s to the quotient state list in the position: quotient_state_position
        quotient_states_list.insert(quotient_states_list.begin() + quotient_state_position, s);
    }
    else{  // we have to insert full_hamiltonian_state_index_list and monomer_state_index_list into our quotient state which is already in the quotient_state_list.
        monomer_state_index_position = binary_insert_dxindex(quotient_states_list[quotient_state_position].monomer_state_index_list, monomer_state_index);

        quotient_states_list[quotient_state_position].full_hamiltonian_state_index_list.insert(quotient_states_list[quotient_state_position].full_hamiltonian_state_index_list.begin() + monomer_state_index_position, full_hamiltonian_state_index );  // quotient_states_list[quotient_state_position] play same role as s.
    }
}

void save_detector_quotient_state_data( const vector <quotient_state> & dlist, ofstream & save){
    int listsize= dlist.size();
    int moddim= dlist[0].vmode.size();
    int xindex_size;
    int q_index_list_size;
    int q_index_size=5;  // size of q_index. (i,j,k,l,p) see construct_q_index() function.
    int i,j,k;
    save<< listsize << " "<< moddim <<endl;
    for(i=0;i<listsize;i++){
        // go through quotient_state d
        // exciton_state_index_list.  vmode.  x_index_length, full_hamiltonian_state_index_list,  monomer_state_index_list.  anharmonic_coupling_info_index_list
        save << dlist[i].exciton_state << " ";
        for(j=0;j<moddim;j++){
            save<< dlist[i].vmode[j]<<" ";
        }
        xindex_size=dlist[i].full_hamiltonian_state_index_list.size();
        save<< xindex_size<<" ";
        for(j=0;j<xindex_size;j++){
            save << dlist[i].full_hamiltonian_state_index_list[j] << " ";
        }
        for(j=0;j<xindex_size;j++){
            save << dlist[i].monomer_state_index_list[j] << " ";
        }
        // output anharmonic_coupling_info_index_list.
        q_index_list_size= dlist[i].anharmonic_coupling_info_index_list.size();
        save<< q_index_list_size<<" ";
        for(j=0;j<q_index_list_size;j++){
            for(k=0;k<q_index_size;k++){
                save << dlist[i].anharmonic_coupling_info_index_list[j][k] << " ";
            }
        }
        save<<endl;
    }
}

void save_detector_quotient_state_data_for2( const vector <quotient_state> & d1list, const vector <quotient_state> & d2list, string path){
    // wrapper function to encode 2 list's data
    ofstream save(path + "quotient_state.txt");
    save_detector_quotient_state_data(d1list,save);
    save_detector_quotient_state_data(d2list,save);
    save.close();
}

void load_detector_quotient_state_data(vector <quotient_state> & dlist, ifstream & load){
    int listsize;
    int moddim;
    string ss;
    int i,j,k;
    int sys_state;
    int xindex_length;
    int q_index_list_size;
    int q_index_size=5;
    int var;
    load>> listsize >> moddim;
    std::getline(load, ss);
    for(i=0;i<listsize;i++) {
        vector<int> vmode;
        load >> sys_state;
        for (j = 0; j < moddim; j++) {
            load >> var;
            vmode.push_back(var);
        }
        quotient_state s(vmode, sys_state);
        load >> xindex_length;
        for (j = 0; j < xindex_length; j++) {
            load >> var;
            s.full_hamiltonian_state_index_list.push_back(var);
        }
        for (j = 0; j < xindex_length; j++) {
            load >> var;
            s.monomer_state_index_list.push_back(var);
        }
        load >> q_index_list_size;
        for (j = 0; j < q_index_list_size; j++) {
            vector<int> q_index;
            for (k = 0; k < q_index_size; k++) {
                load >> var;
                q_index.push_back(var);
            }
            s.anharmonic_coupling_info_index_list.push_back(q_index);
        }
        dlist.push_back(s);
        std::getline(load, ss);
    }
}

void load_detector_quotient_state_data_for2( vector<quotient_state> & d1list, vector <quotient_state> &d2list, string path){
    ifstream load (path + "quotient_state.txt");
    load_detector_quotient_state_data(d1list,load);
    load_detector_quotient_state_data(d2list,load);
    load.close();
}
