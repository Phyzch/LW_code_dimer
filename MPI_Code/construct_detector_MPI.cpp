//
// Created by phyzch on 6/23/20.
// This file contain function that use our previously coded programm to construct matrix.
//
#include"../system.h"
#include"../util.h"
using namespace std;
void Broadcast_dmat_vmode(int stlnum, vector<double> & dmat0,  vector<double> & dmat1,  vector<vector<int>> & vmode0, vector<vector<int>> & vmode1);

// initialize parameters for detector and allocate space there.
void detector::allocate_space() {
    int i, j;

    nmodes = new int[electronic_state_num];  //  number of modes in each detector

    proptime = new double[electronic_state_num];  //  pre-coupling propagation time

    nmax = new int *[electronic_state_num];  // maximum number of quanta in eqch mode.

    dmatsize = new int[electronic_state_num];   // size of detector matrix

    mfreq = new double *[electronic_state_num];  // mfreq: frequency of each mode here.

    electron_phonon_coupling = new double * [electronic_state_num];

    aij = new double * [electronic_state_num];

    dmat = new vector <double> [electronic_state_num];

    dirow = new vector<int> [electronic_state_num];

    dicol = new vector<int> [electronic_state_num];

    dv = new vector<vector<int>> [electronic_state_num];

    // matrix element number for detector matrix
    dmatnum = new int[electronic_state_num];
    // off diagonal matrix element number for detector matrix
    doffnum = new int[electronic_state_num];

    xd_all = new double * [electronic_state_num];
    yd_all = new double * [electronic_state_num];

    // tell detector total matrix size..
    total_dmat_size.reserve(2);
    total_dmat_num.reserve(2);
    total_dmat_off_num.reserve(2);
    dmatsize_each_process = new int * [2];
    dmatsize_offset_each_process = new int * [2];
    doffnum_each_process= new int * [2];
    dmatnum_each_process= new int * [2];;  // record detector matrix element number in each process.
    dmat_offset_each_process= new int * [2];; // record local first detector matrix's index in global matrix.
    for(i=0;i<electronic_state_num;i++){
        dmatsize_each_process[i]= new int [num_proc];
        dmatsize_offset_each_process[i] = new int [num_proc];
        doffnum_each_process[i]= new int [num_proc];
        dmatnum_each_process[i] = new int [num_proc];
        dmat_offset_each_process[i] = new int [num_proc];
    }

};

void detector:: allocate_space_single_detector(int detector_index){
    int i= detector_index;
    nmax[i] = new int[nmodes[i]];
    mfreq[i] = new double[nmodes[i]];
    electron_phonon_coupling[i] = new double [nmodes[i]];
    aij[i] = new double[nmodes[i]];
}

// read parameters for detector matrix construction.
void detector::read_MPI(ifstream & input, ofstream & output, int electronic_state_num1, string path) {
    int i, j;
    double geometric_mean_frequency;

    electronic_state_num = electronic_state_num1;

    allocate_space();
    if (my_id==0){
        // maxdis: maximum allowed distance.
        // cutoff: cutoff strength for intra-detector coupling strength.  cutoff2: cutoff strength for inter-detector coupling.
        input >>  maxdis >> cutoff >> Franck_condon_factor_cutoff;
        output << "Detector  " << maxdis << " " << cutoff <<  endl;
        for (i = 0; i < electronic_state_num; i++) {
            input >> nmodes[i] >> proptime[i];
            allocate_space_single_detector(i);
            for (j = 0; j < nmodes[i]; j++) {
                input >> mfreq[i][j] >> nmax[i][j]  >> electron_phonon_coupling[i][j] ;
            }
        }


        for(i=0; i < electronic_state_num; i++){
            // geometric mean frequency
            geometric_mean_frequency = 1;
            for(j=0; j<nmodes[i];j++){
                geometric_mean_frequency = geometric_mean_frequency * pow(double(mfreq[i][j]) , 1 / double(nmodes[i]));
            }


            output << nmodes[i] << " " << proptime[i] << endl;
            for(j=0;j<nmodes[i];j++){
                aij[i][j] = pow(mfreq[i][j] , 0.5)/ 270;
                // aij corresponding to scaling factor for f= f_{bright}/2 cm^{-1}.
                output << mfreq[i][j] << " " << nmax[i][j] << endl;
            }
        }
    }

    MPI_Bcast(&nmodes[0], electronic_state_num, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&proptime[0], electronic_state_num, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if(my_id!=0){
        for(i=0; i < electronic_state_num; i++){
            allocate_space_single_detector(i);
        }
    }

    // temporary value for calculation of coupling.
    deln = new int[max(nmodes[0],nmodes[1])];
    nbar = new double[max(nmodes[0],nmodes[1])];

    MPI_Bcast(&maxdis,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&cutoff,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&Franck_condon_factor_cutoff, 1, MPI_DOUBLE, 0,  MPI_COMM_WORLD);

    for(i=0; i < electronic_state_num; i++){
        MPI_Bcast(&mfreq[i][0],nmodes[i],MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Bcast(&electron_phonon_coupling[i][0] ,nmodes[i],MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Bcast(&nmax[i][0],nmodes[i],MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(&aij[i][0],nmodes[i],MPI_DOUBLE,0,MPI_COMM_WORLD);
    }
};

void detector:: construct_dmatrix_MPI(ifstream & input, ofstream & output, ofstream & log, vector<double> & dmat_diagonal_global0, vector<double> & dmat_diagonal_global1,  vector<vector<int>> & vmode0, vector<vector<int>> & vmode1) {
    int m, i;

    construct_dv_dirow_dicol_dmatrix_MPI(log, dmat_diagonal_global0, dmat_diagonal_global1, vmode0, vmode1);

    compute_important_state_index();
    // -------------------------- constructing off-diagonal term for detector  -----------------------------
    compute_detector_offdiag_part_MPI(log, dmat_diagonal_global0, dmat_diagonal_global1, vmode0, vmode1);

    broadcast_dmatnum_doffnum();
    broadcast_total_dmat();

    if(my_id==0){
        for(m=0; m < electronic_state_num; m++){
            output << "Matrix for detector " << m << " : " << total_dmat_size[m] << "  " << total_dmat_off_num[m] << endl;
        }

    }

    xd = new vector <double> [electronic_state_num];
    yd = new vector<double> [electronic_state_num];
    for (i = 0; i < electronic_state_num; i++) {
        xd[i].reserve(dmatsize[i]);
        yd[i].reserve(dmatsize[i]);
    }
}

void detector:: construct_dv_dirow_dicol_dmatrix_MPI(ofstream & log,vector<double> & dmat0,  vector<double> & dmat1,  vector<vector<int>> & vmode0, vector<vector<int>> & vmode1){
    int i,m;
    int vsize,vsize2;
    int displacement;
    if(my_id==0){
        total_dmat_size[0]= dmat0.size();
        if(electronic_state_num == 2){
            total_dmat_size[1] = dmat1.size();
        }
        else{
            total_dmat_size[1]=0;
        }
    }
    MPI_Bcast(&total_dmat_size[0],2,MPI_INT,0,MPI_COMM_WORLD);
    for(i=0; i < electronic_state_num; i++){
        xd_all[i] = new double [total_dmat_size[i]];
        yd_all[i] = new double [total_dmat_size[i]];
    }
    vector <int> * dirow_all, *dicol_all;
    vector< double> * dmat_all;
    int ** vector_size, ** displacement_list;

    dv_all = new vector<vector<int>>[2];
    dmat_all = new vector<double>[2];  // do not confuse dmat_all with total_dmat. dmat_all only contain diagonal term
    dirow_all = new vector<int>[2];
    dicol_all = new vector<int>[2];
    vector_size = new int *[2];
    displacement_list = new int *[2];

    Broadcast_dmat_vmode(electronic_state_num, dmat0, dmat1, vmode0, vmode1); // broadcast dmat_diagonal_global0, dmat_diagonal_global1, vmode0, vmode1 to all process to compute off-diagonal matrix.
    dv_all[0] = vmode0;
    dv_all[1] = vmode1;

    if(my_id==0) {
        dmat_all[0] = dmat0;
        dmat_all[1] = dmat1;
        // prepare dirow, dicol for broadcasting.
        for (m = 0; m < electronic_state_num; m++) {
            int size = dmat_all[m].size();
            for (i = 0; i < size; i++) {
                dirow_all[m].push_back(i);
                dicol_all[m].push_back(i);
            }
        }
        // prepare vector size and vector displacement:
        for (m = 0; m < electronic_state_num; m++) {
            vsize = total_dmat_size[m] / num_proc;
            vsize2 = total_dmat_size[m] - (num_proc - 1) * vsize;
            vector_size[m] = new int[num_proc];  // size of vector to scatter to each process.
            displacement_list[m] = new int[num_proc]; // displacement of vector to scatter to each process.
            displacement = 0;
            for (i = 0; i < num_proc - 1; i++) {
                vector_size[m][i] = vsize;
                displacement_list[m][i] = displacement;
                displacement = displacement + vsize;
            }
            vector_size[m][num_proc - 1] = vsize2;
            displacement_list[m][num_proc - 1] = displacement;
        }
    }
    Scatter_dirow_dicol_dmatrix(dmat_all, dirow_all, dicol_all, vector_size, displacement_list, log);
    Scatter_dv(total_dmat_size);  // scatter dv_all
    // construct detector matrix size for each process.
    if(my_id != num_proc-1){
        for(m=0; m < electronic_state_num; m++) {
            vsize= total_dmat_size[m]/num_proc;
            dmatsize[m] = vsize;
        }
    }
    else{
        for(m=0; m < electronic_state_num; m++){
            vsize= total_dmat_size[m] - total_dmat_size[m]/num_proc *(num_proc-1);
            dmatsize[m]=vsize;
        }
    }
    for(m=0; m < electronic_state_num; m++){
        MPI_Allgather(&dmatsize[m],1, MPI_INT,&dmatsize_each_process[m][0],1,MPI_INT,MPI_COMM_WORLD);
    }

    for(m=0;m<electronic_state_num;m++){
        dmatsize_offset_each_process[m][0] = 0;
        for(i=1;i<num_proc;i++){
            dmatsize_offset_each_process[m][i] = dmatsize_offset_each_process[m][i-1] + dmatsize_each_process[m][i-1];
        }
    }

    if(my_id == 0) {
        for (m = 0; m < electronic_state_num; m++) {
            delete[] vector_size[m];
            delete[] displacement_list[m];
        }
    }
    delete[] vector_size;
    delete[] displacement_list;
    delete[] dirow_all;
    delete[] dicol_all;
}

void Broadcast_dmat_vmode(int stlnum, vector<double> & dmat0,  vector<double> & dmat1,  vector<vector<int>> & vmode0, vector<vector<int>> & vmode1){
    int i,j,m;
    int dmat0_size, dmat1_size;
    // we broadcast dmat_diagonal_global0, dmat_diagonal_global1, .. to all other process. This is need for us to compute off diagonal matrix
    // first allocate space for dmat_diagonal_global0 , dmat_diagonal_global1.
    if(my_id==0){
        dmat0_size= dmat0.size();
        dmat1_size= dmat1.size();
    }
    MPI_Bcast(&dmat0_size,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&dmat1_size,1,MPI_INT,0,MPI_COMM_WORLD);
    if(my_id!=0) {
        dmat0.resize(dmat0_size);
        dmat1.resize(dmat1_size);
    }
    // Broadcast dmat_diagonal_global0, dmat_diagonal_global1
    MPI_Bcast((void *) dmat0.data(),dmat0_size,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast((void *) dmat1.data(),dmat1_size,MPI_DOUBLE,0,MPI_COMM_WORLD);

    // ----------------Broadcast vmode0, vmode1 -------------------------------
    vector<vector<int>> * v_ptr;
    // ------------- For sending vmode0,vmode1 --------------------------
    vector<int> vmode_1d;
    vector<int> displacement;
    vector <int> element_size;
    int vmode_1d_size;
    int element_number;
//-------------------------------------------
    for (m = 0; m < stlnum; m++) {
        vmode_1d.clear();
        displacement.clear();
        element_size.clear();
        if (m == 0) v_ptr = &(vmode0);
        else v_ptr = &(vmode1);
        if(my_id==0) {
            convert_dv(*v_ptr, vmode_1d, displacement, element_size);
            vmode_1d_size = vmode_1d.size();
            MPI_Bcast(&vmode_1d_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
            element_number = element_size.size();
            MPI_Bcast(&element_number, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast((void *) element_size.data(), element_number, MPI_INT, 0,
                      MPI_COMM_WORLD); // send size of each 2d element
            MPI_Bcast((void *) vmode_1d.data(), vmode_1d_size, MPI_INT, 0,
                      MPI_COMM_WORLD); // send converted 1d data to other process.
        }
        else{
            MPI_Bcast(&vmode_1d_size,1,MPI_INT,0,MPI_COMM_WORLD);
            vmode_1d.reserve(vmode_1d_size);
            MPI_Bcast(&element_number,1,MPI_INT,0,MPI_COMM_WORLD);
            element_size.reserve(element_number);
            MPI_Bcast( (void *) element_size.data(), element_number,MPI_INT,0,MPI_COMM_WORLD);
            MPI_Bcast((void *) vmode_1d.data(),vmode_1d_size,MPI_INT,0,MPI_COMM_WORLD);

            int index=0;
            // convert vmode_1d to vmode0 / vmode1
            for(i=0;i<element_number;i++){
                vector<int> dv_vmode;
                dv_vmode.reserve(element_size[i]);
                for(j=0;j<element_size[i];j++){
                    dv_vmode.push_back(vmode_1d[index]);
                    index++;
                }
                (*v_ptr).push_back(dv_vmode);
            }
        }
    }
}




void detector::compute_detector_offdiag_part_MPI(ofstream & log,vector<double> & dmat0,  vector<double> & dmat1,  vector<vector<int>> & vmode0, vector<vector<int>> & vmode1){
    int i,j,m,k;
    int begin_index;
    int ntot;
    double value, lij;
    vector<vector <int>> * vmode_ptr;
    vector<double> * dmat_ptr;
    bool exist;
    int position;
    double random_number;
    // different process do different amount of work.
    for(m=0; m < electronic_state_num; m++){
        if(m==0){
            vmode_ptr = &(vmode0);
            dmat_ptr= &(dmat0);
        }
        else {
            vmode_ptr= &(vmode1);
            dmat_ptr= &(dmat1);
        }
        begin_index= total_dmat_size[m]/num_proc * my_id;
        // compute off diagonal matrix element, using dmat_diagonal_global0, dmat_diagonal_global1.
        for(i=begin_index;i<begin_index + dmatsize[m];i++){  // for my_id==0 , O(dmatsize * dmatsize/ proc_num)
            for(j=0;j<total_dmat_size[m];j++){ // j is different from serial program. Now we record both symmetric Hamiltonian element
                if (i==j) continue;
                ntot=0;
                for(k=0;k<nmodes[m];k++){
                    deln[k]= abs( (*vmode_ptr)[i][k] - (*vmode_ptr)[j][k] ); // same as deln[k] = abs(dv[m][i][k] - dv[m][j][k]);
                    nbar[k]= sqrt(sqrt(double(max(1, (*vmode_ptr)[i][k])) * double(max(1, (*vmode_ptr)[j][k]  ))));
                    ntot= ntot+ deln[k];
                }
                // this is because we don't have 1st order ladder operator in Harmonic oscillator's expression
                if (ntot == 2) {
                    for (k = 0; k < nmodes[m]; k++) {
                        if (deln[k] == 2) deln[k] = 4;
                        if (deln[k] == 1) deln[k] = 2;
                    }
                } else if (ntot == 1) {
                    for (k = 0; k < nmodes[m]; k++) {
                        if (deln[k] == 1) deln[k] = 3;
                    }
                } else if (ntot == 0) {
                    log << "Error! The off-diagonal element must differ in q.n." << endl;
                    MPI_Abort(MPI_COMM_WORLD,-8);
                }
                // check ntot with maxdis in quantum number space:
                if (ntot < maxdis) {

                    if (ntot % 2 == 0) {
                        value = V_intra;  // V=0.03 as requirement.
                    } else {
                        value = -V_intra;
                    }
                    for (k = 0; k < nmodes[m]; k++) {
                        value = value * pow(aij[m][k]* nbar[k], deln[k]);
                    }
                    if ( (*dmat_ptr)[i] != (*dmat_ptr)[j] ) {
                        lij = abs(value / ((*dmat_ptr)[i] - (*dmat_ptr)[j]));
                        if (lij > cutoff) {
                            dmat[m].push_back(value);
                            dirow[m].push_back(i);
                            dicol[m].push_back(j);
                        }
                    } else {
                        dmat[m].push_back(value);
                        dirow[m].push_back(i);
                        dicol[m].push_back(j);
                    }
                }
            }
        }
    }
}

void detector:: broadcast_dmatnum_doffnum(){
    int m,i;
    // compute dmatnum and doffnum and dmatnum_each_process, total_dmatnum, total_doffnum etc.
    for(m=0; m < electronic_state_num; m++){
        dmatnum[m]= dmat[m].size();
        doffnum[m]=dmatnum[m] - dmatsize[m];
    }
    // compute toatl_dmatnum, total_dmatoffnum
    for(m=0; m < electronic_state_num; m++){
        MPI_Allreduce(&dmatnum[m],&total_dmat_num[m],1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
        MPI_Allreduce(&doffnum[m],&total_dmat_off_num[m],1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    }
    // broadcast detector matrix number and detector matrix off-diagonal number in each process.
    for(m=0; m < electronic_state_num; m++){
        MPI_Allgather(&dmatnum[m],1,MPI_INT,&dmatnum_each_process[m][0],1,MPI_INT,MPI_COMM_WORLD);
        MPI_Allgather(&doffnum[m], 1, MPI_INT, &doffnum_each_process[m][0], 1, MPI_INT,MPI_COMM_WORLD);
    }
    // compute offset of detector matrix for each process.
    for(m=0; m < electronic_state_num; m++){
        dmat_offset_each_process[m][0]=0;
        for(i=1;i<num_proc;i++){
            dmat_offset_each_process[m][i] = dmat_offset_each_process[m][i-1] + dmatnum_each_process[m][i-1];
        }
    }
}


void detector:: broadcast_total_dmat(){
    /*
     construct total_dmat, total_dirow, total_dicol
     use dmatnum_each_process,  dmat_offset_each_process\
    */
    total_dmat= new double * [electronic_state_num];
    total_dirow= new int * [electronic_state_num];
    total_dicol= new int * [electronic_state_num];
    int m;
    for(m=0; m < electronic_state_num; m++){
        total_dmat[m] = new double [total_dmat_num[m]];
        total_dirow[m] = new int [total_dmat_num[m]];
        total_dicol[m] = new int [total_dmat_num[m]];
        MPI_Allgatherv(&dmat[m][0],dmatnum[m],MPI_DOUBLE,
                &total_dmat[m][0],dmatnum_each_process[m],dmat_offset_each_process[m],MPI_DOUBLE,MPI_COMM_WORLD);
        MPI_Allgatherv(&dirow[m][0],dmatnum[m],MPI_INT,
                &total_dirow[m][0],dmatnum_each_process[m],dmat_offset_each_process[m],MPI_INT,MPI_COMM_WORLD);
        MPI_Allgatherv(&dicol[m][0],dmatnum[m],MPI_INT,
                &total_dicol[m][0],dmatnum_each_process[m],dmat_offset_each_process[m],MPI_INT,MPI_COMM_WORLD);
    }
}

//-------------------------------- Construct detector wave function ------------------------------
void detector::initialize_detector_state_MPI(ofstream & log) {
    // this is the code for initializing detector state's wavefunction.
    int m, i,j;
    bool check_mark;
    double norm;
    double total_norm;
    // this way we initialize wave function according to distribution
    for (m = 0; m < electronic_state_num; m++) {
        norm = 0;
        int ** control_state;
        control_state=initial_detector_state;
        for (i = 0; i < dmatsize[m]; i++) {
            check_mark = true;  // check if state is initial_detector state/ bright state.

            for(j=0;j<nmodes[m];j++){
                if(dv[m][i][j] != control_state[m][j]){
                    check_mark=false;
                    break;
                }
            }

            if(check_mark) {
                xd[m].push_back(1);
                yd[m].push_back(0);
                norm = norm + pow(xd[m][i], 2) + pow(yd[m][i],2);
            }

            else{
                xd[m].push_back(0);
                yd[m].push_back(0);
            }

        }
        MPI_Allreduce(&norm,&total_norm,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        //  ---------------------------------------------------------------
        if(total_norm==0){
            if(my_id==0) {
                cout << "Norm for detector state is 0" << endl;
                log << "Norm for detector state is 0" << endl;
                MPI_Abort(MPI_COMM_WORLD,-10);
            }
        }

        total_norm=1/sqrt(total_norm);
        for (i = 0; i < dmatsize[m]; i++) {
            xd[m][i] = xd[m][i] * total_norm;
        }
    }


}

void detector::construct_initial_state_MPI( vector<vector<int>> & initial_state_quantum_number){
    // MPI version of construct initial state
    /*
     *  Initial detector state: detector state populated at beginning of simulation
     *  Initial_Detector_energy: Energy of detector at initial state.
     */
    int m,i;

    initial_detector_state= new int * [electronic_state_num];
    initial_state_energy= new double [electronic_state_num];
    for(m=0; m < electronic_state_num; m++){
        initial_detector_state[m]= new int [nmodes[m]];
        initial_state_energy[m]=0;
    }
    if(my_id==0){
        for(m=0; m < electronic_state_num; m++) {
            for (i = 0; i < nmodes[m]; i++) {
                initial_detector_state[m][i] = initial_state_quantum_number[m][i];
            }
        }
    }
    for(m=0; m < electronic_state_num; m++){ // Broad cast initial detector state.
        MPI_Bcast(&initial_detector_state[m][0],nmodes[m],MPI_INT,0,MPI_COMM_WORLD);
    }

    for(m=0; m < electronic_state_num; m++){
        // initialize our initial detector state.  set dark mode's quanta equal to bright_state.
        for(i=0;i<nmodes[m];i++){
            initial_state_energy[m]= initial_state_energy[m] + initial_detector_state[m][i] * mfreq[m][i] - pow( mfreq[m][i] * initial_detector_state[m][i] , 2)/ (4 * self_anharmonicity_D);
        }
    }

    if(my_id==0){  // output initial detector state to output.txt
        cout <<"Initial detector state:"<<endl;
        for(m=0; m < electronic_state_num; m++){
            for(i=0;i<nmodes[m];i++){
                cout<< initial_detector_state[m][i] <<" ";
            }
            cout<<endl;
        }
    }
}

void detector:: compute_important_state_index(){
    // compute bright state index and initial state index for two detector.
    int m,i;

    vector<int> special_state_vec ;
    int ** special_state;
    int * special_state_index;
    int * special_state_pc_id;
    int position;
    bool exist;
    bool * exist_bool_for_pc = new bool [num_proc];
    initial_state_index = new int [2];
    initial_state_pc_id = new int [2];
    // we initialize bright_state_index and initial_state_index.
    // We do this because sometimes we may not include bright state in our simulation, then this index need to be initialized.
    for(m=0; m < electronic_state_num; m++){
        initial_state_index[m]=0;
        initial_state_pc_id[m] = 0;
    }

    // loop for bright state and initial state
    special_state = initial_detector_state;
    special_state_index = initial_state_index;
    special_state_pc_id = initial_state_pc_id;
    // loop for two detector.
    for (m = 0; m < electronic_state_num; m++) {
        special_state_vec.clear();
        for (i = 0; i < nmodes[m]; i++) {
            special_state_vec.push_back(special_state[m][i]);
        }
        position=find_position_for_insert_binary(dv[m],special_state_vec,exist);
        MPI_Allgather(&exist,1,MPI_C_BOOL,&exist_bool_for_pc[0],1,MPI_C_BOOL,MPI_COMM_WORLD);
        special_state_pc_id [m] = -1;
        for(i=0;i<num_proc;i++){
            if(exist_bool_for_pc[i]){
                special_state_pc_id[m] = i;
            }
        }
        if(special_state_pc_id[m] == -1){
            if(my_id==0){
                cout<<" Can not find initial state or brigth state in all vmode. Must have bug here."<<endl;
                MPI_Abort(MPI_COMM_WORLD,-7);
            }
        }
        MPI_Bcast(&position,1, MPI_INT,special_state_pc_id[m],MPI_COMM_WORLD);
        special_state_index [m] = position;
    }

    delete [] exist_bool_for_pc;
}

