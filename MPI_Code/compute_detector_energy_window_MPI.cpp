//
// Created by phyzch on 7/22/20.
//
# include"../util.h"
# include "../system.h"

using namespace std;
int Rmax; // maximum distance allowed in monomer state space.


int compute_state_space_distance(const vector<int> & ndetector, int * state1, int moddim){
    // compute distance between two state : ndetector and state1. Dimension of mod is moddim
    int i;
    int distance=0;
    for(i=0;i<moddim;i++){
        distance= distance + abs(state1[i] - ndetector[i]);
    }
    return distance;
}


int vmode_compare(const vector <int> & lhs, const vector <int> & rhs){
    // compare function used for find_position_for_insert_binary.
    int i,j;
    int size= lhs.size();
    if(size != rhs.size()){
        cout<<"size of vector do not match."<<endl;
        exit(-5);
    }
    for(i=0;i<size;i++){
        if(lhs[i] > rhs[i]) return 1;
        else if(lhs[i] < rhs[i]) return -1;
    }
    return 0;
}

//called in compute_matrix_size() function.
void copy_array(vector<int> & temporary_vector, const vector<int> & ndetector,int modnumber) {
    // code used in compute_matrix_size() function
    // modnumber here is the size of ndetector (the two detectors have the same size)
    int i;
    for (i = 0; i < modnumber; i++) {
        temporary_vector[i] = ndetector[i];
    }
}


// find position to insert the corresponding monomer state, called in compute_matrix_size() function.  Binary search.
int find_position_for_insert_binary(const vector<vector<int>> & vmode, const vector<int> & ndetector, bool & exist) {
    // code used in compute_matrix_size() function
    // first compare first mode, 0 is in front of 1. So finally our mode is |00000> , |00001>,|00010>, etc, |10000>
    if (vmode.size() ==0){
        exist = false;
        return 0;
    }
    int left_flag=0;
    int right_flag=vmode.size() -1;
    int position = (left_flag + right_flag) /2;
    exist = false;
    int mark;
    // binary search using compare function
    while(right_flag>left_flag) {
        mark = vmode_compare(ndetector, vmode[position]);
        if (mark == 1) {
            // ndetector should be after position
            left_flag = position+1;
            position = (left_flag + right_flag) / 2;
        }
        else if ( mark == -1) {
            right_flag = position-1;
            position = (left_flag + right_flag) / 2;
        }
        else if (mark== 0) {
            exist = true;
            return position;
        }
    }
    // now left == right. Decide position to insert and return it.
    mark = vmode_compare(ndetector, vmode[position]);
    if( mark == 0){
        exist=true;
        return position;  // we will not insert.
    }
    else if ( mark == -1){
        exist=false;
        return position;   // we will insert before this position.
    }
    else if ( mark == 1){
        exist=false;
        return position +1; // we will insert after this position.
    }

    //
}


void full_system:: compute_monomer_vib_state_basis_set_size_MPI( ){
    if(my_id==0){
        int i, j, k;
        int i1, i2;
        double  monomer0_energy, monomer1_energy;
        // monomer0_qn and monomer1_qn indicate current monomer mode index we are calculating energy.
        vector<int> monomer0_qn(d.nmodes[0]);
        vector <int> monomer1_qn(d.nmodes[1]);

        // record size of total matrix
        int location;
        bool exist=false;

        int state_space_distance;

        // --------------- code for constructing basis set for the first monomer ---------------------------
        monomer0_qn[0] = -1; // this is for:  when we go into code: monomer0_qn[i]= monomer0_qn[i]+1, our first state is |000000>
        while (1) {
            label2:;  // label2 is for detector1 to jump out of while(1) loop (this is inner layer of while(1))
            monomer0_energy = 0;
            for (i1 = 0; i1 < d.nmodes[0]; i1++) {  // loop through detector0
                // define the way we loop through detector0:
                monomer0_qn[i1] = monomer0_qn[i1] + 1;
                if (monomer0_qn[i1] <= d.nmax[0][i1]) break;
                if (monomer0_qn[d.nmodes[0] - 1] > d.nmax[0][d.nmodes[0] - 1]) {
                    monomer0_qn[d.nmodes[0] - 1] = 0;
                    goto label1;  // use goto to jump out of nested loop
                }
                monomer0_qn[i1] = 0;
            }
            // calculate monomer 0 energy
            for (i = 0; i < d.nmodes[0]; i++) {
                if (self_anharmonicity_bool){
                    // add self-anharmonicity
                    monomer0_energy = monomer0_energy + d.mfreq[0][i] * (monomer0_qn[i] - pow(monomer0_qn[i], 2) * d.mfreq[0][i] / (4 * self_anharmonicity_D) );
                }
                else{
                    monomer0_energy = monomer0_energy + d.mfreq[0][i] * monomer0_qn[i] ;
                }
            }

            //--------------------------------------------------------------------------------------------
            // criteria below make sure monomer 0 's energy is reasonable.
            if (monomer0_energy > d.initial_state_energy[0] + d.vibrational_energy_window_size) {
                // monomer 0's energy can not be larger than its initial energy + photon energy
                // jump to next monomer state.
                k = 0;
                while (monomer0_qn[k] == 0) {
                    monomer0_qn[k] = d.nmax[0][k];
                    k++;
                    if (k >= d.nmodes[0]) {
                        break;
                    }
                }
                if (k < d.nmodes[0]) {
                    monomer0_qn[k] = d.nmax[0][k];
                }
                goto label2;
            }

            // criteria for energy window around bright_state and lower bright state for monomer 0
            if ((monomer0_energy > d.initial_state_energy[0] - d.vibrational_energy_window_size and
                 monomer0_energy < d.initial_state_energy[0] + d.vibrational_energy_window_size)
                    )
                ;
            else {
                goto label2;
            }

            //------------------------------------------------------------------------------------------------
            // criteria below make sure monomer 1 can not be too far away from bright state and lower bright state.
            state_space_distance =
                    compute_state_space_distance(monomer0_qn, d.initial_vibrational_state[0], d.nmodes[0]);

            // we use distance constraint for state whose energy is between two
            if (state_space_distance > Rmax) {
                goto label2;
            }

            //--------------------------------------insert this state in monomer's state.-----------------------------------------------------------
            location = find_position_for_insert_binary(monomer_qn_list0, monomer0_qn, exist);  // we check if this mode exist and the location we have to insert this state at the same time.
            if (!exist) {
                if(location > monomer_qn_list0.size()){
                    cout << "Something wrong when inserting element into monomer_qn_list0  "<<endl;
                    MPI_Abort(MPI_COMM_WORLD,-15);
                }
                // when we push back we should consider arrange them in order. We compute location to insert in find_position_for_insert_binary() function:
                monomer_qn_list0.insert(monomer_qn_list0.begin() + location, monomer0_qn);
                monomer1_vib_state_energy_all_pc.insert(monomer1_vib_state_energy_all_pc.begin() + location, monomer0_energy);
            }
        }
        label1:;



        // ---------- repeat the code for second monomer -----------------------------
        monomer1_qn[0] = -1; // this is when we go into code: monomer1_qn[i] = monomer1_qn[i]+1. our first state is |000000>
        while (1) { // loop through monomer 1
            label3:;
            monomer1_energy = 0;
            for (i2 = 0; i2 < d.nmodes[1]; i2++) {
                // define the way we loop through detector1
                monomer1_qn[i2] = monomer1_qn[i2] + 1;
                if (monomer1_qn[i2] <= d.nmax[1][i2]) break;
                if (monomer1_qn[d.nmodes[1] - 1] > d.nmax[1][d.nmodes[1] - 1]) {
                    monomer1_qn[d.nmodes[1] - 1] = 0;
                    goto label4;
                }
                monomer1_qn[i2] = 0;
            }
            // calculate monomer 1 energy
            for (i = 0; i < d.nmodes[1]; i++) {
                if (self_anharmonicity_bool){
                    // add self-anharmonicity
                    monomer1_energy = monomer1_energy + d.mfreq[1][i] * (monomer1_qn[i] - pow(monomer1_qn[i], 2) * d.mfreq[1][i] / (4 * self_anharmonicity_D) );
                }
                else{
                    monomer1_energy = monomer1_energy + monomer1_qn[i] * d.mfreq[1][i];
                }
            }
            // --------------------------------------------------------------
            //  criteria below make sure monomer 1's energy is reasonable:
            //  we exclude states whose vibrational energy is too high.
            if (monomer1_energy > d.initial_state_energy[1] + d.vibrational_energy_window_size) {
                // initial energy is system energy.
                // monomer 1 's energy can not be larger than its initial energy + photon energy
                j = 0;
                while (monomer1_qn[j] == 0) { // go to first mode whose n!=0;
                    monomer1_qn[j] = d.nmax[1][j];
                    j++;
                    if (j >= d.nmodes[1]) {
                        break;
                    }
                }
                if (j < d.nmodes[1]) {
                    monomer1_qn[j] = d.nmax[1][j];
                }
                goto label3;
            }

            // criteria for energy window around bright_state and lower bright state for monomer 1
            if ((monomer1_energy > d.initial_state_energy[1] -
                                   d.vibrational_energy_window_size and
                 monomer1_energy < d.initial_state_energy[1] +
                                   d.vibrational_energy_window_size )
                    )
            {  // criteria here means we only consider monomer state whose energy is within small energy window
                ;
            } else {
                goto label3;
            }
            //------------------------------------------------------------------------------------------------
            // criteria below make sure monomer 1 can not be too far away from bright state and lower bright state.

            state_space_distance = compute_state_space_distance(monomer1_qn, d.initial_vibrational_state[1],
                                                                d.nmodes[1]);

            if (state_space_distance > Rmax ) {
                goto label3;
            }
            location = find_position_for_insert_binary(monomer_qn_list1, monomer1_qn, exist);
            if (!exist) {
                if(location > monomer_qn_list1.size()){
                    cout <<"Something wrong when inserting into monomer_qn_list1. " <<endl;
                    MPI_Abort(MPI_COMM_WORLD,-15);
                }
                monomer_qn_list1.insert(monomer_qn_list1.begin() + location, monomer1_qn);
                monomer2_vib_state_energy_all_pc.insert(monomer2_vib_state_energy_all_pc.begin() + location, monomer1_energy);
            }
        }
        label4:;

    }
}