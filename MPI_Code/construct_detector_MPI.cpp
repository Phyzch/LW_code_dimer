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
    dmatdim = detdim*detdim / fillfrac;

    nmodes = new int[electronic_state_num];  //  number of modes in each detector

    proptime = new double[electronic_state_num];  //  pre-coupling propagation time

    nmax = new int *[electronic_state_num];  // maximum number of quanta in eqch mode.

    dmatsize = new int[electronic_state_num];   // size of detector matrix

    mfreq = new double *[electronic_state_num];  // mfreq: frequency of each mode here.

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
    doffnum_each_process= new int * [2];
    dmatnum_each_process= new int * [2];;  // record detector matrix element number in each process.
    dmat_offset_each_process= new int * [2];; // record local first detector matrix's index in global matrix.
    for(i=0;i<electronic_state_num;i++){
        dmatsize_each_process[i]= new int [num_proc];
        doffnum_each_process[i]= new int [num_proc];
        dmatnum_each_process[i] = new int [num_proc];
        dmat_offset_each_process[i] = new int [num_proc];
    }

};

void detector:: allocate_space_single_detector(int detector_index){
    int i= detector_index;
    nmax[i] = new int[nmodes[i]];
    mfreq[i] = new double[nmodes[i]];
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
        input >>  maxdis >> cutoff;
        output << "Detector  " << maxdis << " " << cutoff <<  endl;
        for (i = 0; i < electronic_state_num; i++) {
            input >> nmodes[i] >> proptime[i];
            allocate_space_single_detector(i);
            for (j = 0; j < nmodes[i]; j++) {
                input >> mfreq[i][j] >> nmax[i][j]  ;
            }
        }

        // geometric mean frequency
        geometric_mean_frequency = 1;
        for(j=0; j<nmodes[i];j++){
            geometric_mean_frequency = geometric_mean_frequency * pow(double(mfreq[i][j]) , 1 / double(nmodes[i]));
        }

        for(i=0; i < electronic_state_num; i++){
            output << nmodes[i] << " " << proptime[i] << endl;
            for(j=0;j<nmodes[i];j++){
                aij[i][j] = a_intra * pow(double(mfreq[i][j]) / geometric_mean_frequency ,0.5);
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

    for(i=0; i < electronic_state_num; i++){
        MPI_Bcast(&mfreq[i][0],nmodes[i],MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Bcast(&nmax[i][0],nmodes[i],MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(&aij[i][0],nmodes[i],MPI_DOUBLE,0,MPI_COMM_WORLD);
    }
};

void detector:: construct_dmatrix_MPI(ifstream & input, ofstream & output, ofstream & log,  vector<double> & dmat0,  vector<double> & dmat1,  vector<vector<int>> & vmode0, vector<vector<int>> & vmode1) {
    int m, i;
    construct_dv_dirow_dicol_dmatrix_MPI(log, dmat0, dmat1, vmode0, vmode1);
    compute_important_state_index();
    // -------------------------- Two different way of constructing off-diagonal term for detector  -----------------------------
    // 1:  traditional way of constructing detector off-diagonal part of matrix
    compute_detector_offdiag_part_MPI(log,dmat0,dmat1,vmode0,vmode1);

    //--------------------------------------------------------------------------------------------------
    update_initial_state_energy();

    output_state_density(dmat0,dmat1);

    broadcast_dmatnum_doffnum();
    broadcast_total_dmat();

    if(my_id==0){
        for(m=0; m < electronic_state_num; m++){
            output << "Matrix for detector " << m << " : " << total_dmat_size[m] << "  " << total_dmat_off_num[m] << endl;
        }

    }

    ofstream local_density_of_state_output(path + "local_density_of_state.txt");
    compute_local_density_of_state(local_density_of_state_output,dmat0,dmat1);
    local_density_of_state_output.close();

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

    Broadcast_dmat_vmode(electronic_state_num, dmat0, dmat1, vmode0, vmode1); // broadcast dmat0, dmat1, vmode0, vmode1 to all process to compute off-diagonal matrix.
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
    // we broadcast dmat0, dmat1, .. to all other process. This is need for us to compute off diagonal matrix
    // first allocate space for dmat0 , dmat1.
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
    // Broadcast dmat0, dmat1
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
        // compute off diagonal matrix element, using dmat0, dmat1.
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
    if(electronic_state_num == 1) {
        if(my_id == 0) {
            xd[1].push_back(1);
            yd[1].push_back(0); // in case electronic_state_num==1
        }
    }
    // this way we initialize wave function according to distribution
    for (m = 0; m < electronic_state_num; m++) {
        norm = 0;
        int ** control_state;
        control_state=initial_detector_state;
        for (i = 0; i < dmatsize[m]; i++) {
            check_mark = true;  // check if state is initial_detector state/ bright state.

            for(j=0;j<nmodes[m];j++){
                if(dv[m][i][j]!=control_state[m][j]){
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

void detector::construct_bright_state_MPI(ifstream & input, ofstream & output){
    // MPI version of construct bright state
    /*
     *  Initial detector state: detector state populated at beginning of simulation
     *  Bright  detector state: detector state with bright mode set to 1 ,all other mode equal to initial detector state
     *  Initial_Detector_energy: Energy of detector at initial state.
     */
    int my_id;
    MPI_Comm_rank(MPI_COMM_WORLD,&my_id);
    int m,i,j;

    initial_detector_state= new int * [electronic_state_num];
    initial_state_energy= new double [electronic_state_num];
    for(m=0; m < electronic_state_num; m++){
        initial_detector_state[m]= new int [nmodes[m]];
        initial_state_energy[m]=0;
    }
    if(my_id==0){
        for(m=0; m < electronic_state_num; m++) {
            for (i = 0; i < nmodes[m]; i++) {
                input >> initial_detector_state[m][i];
            }
        }
    }
    for(m=0; m < electronic_state_num; m++){ // Broad cast initial detector state.
        MPI_Bcast(&initial_detector_state[m][0],nmodes[m],MPI_INT,0,MPI_COMM_WORLD);
    }

    for(m=0; m < electronic_state_num; m++){
        // initialize our initial detector state.  set dark mode's quanta equal to bright_state.
        for(i=0;i<nmodes[m];i++){
            initial_state_energy[m]= initial_state_energy[m] + initial_detector_state[m][i] * mfreq[m][i];
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

void detector:: update_initial_state_energy(){
    // we update energy of initial and bright state of detector. since in Van Vleck transformation, the energy level is shhifted.
    int m,i;
    for(m=0; m < electronic_state_num; m++){
        initial_state_energy[m] = 0;
        if(my_id == initial_state_pc_id[m]) {
            initial_state_energy[m] = dmat[m][initial_state_index[m]];
        }
        MPI_Bcast(&initial_state_energy[m], 1, MPI_DOUBLE, initial_state_pc_id[m], MPI_COMM_WORLD);
    }
}
void detector:: compute_important_state_index(){
    // compute bright state index and initial state index for two detector.
    int m,i,j;
    int index;
    bool check_mark1, check_mark2;
    vector<int> special_state_vec ;
    int ** special_state;
    int * special_state_index;
    int * special_state_pc_id;
    int position;
    bool exist;
    bool * exist_bool_for_pc = new bool [num_proc];
    bright_state_index = new int [2];
    initial_state_index = new int [2];
    bright_state_pc_id = new int [2];
    initial_state_pc_id = new int [2];
    // we initialize bright_state_index and initial_state_index.
    // We do this because sometimes we may not include bright state in our simulation, then this index need to be initialized.
    for(m=0; m < electronic_state_num; m++){
        bright_state_index[m] = 0;
        initial_state_index[m]=0;
        bright_state_pc_id[m] = 0;
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


void compute_detector_state_density(vector <int> & state_density, vector <double> & energy_distribution,
                                    vector<double> & energy_range){
    /*
     input: state_density: D_{l}(E) we are going to calculate
            energy_distribution: energy of different state in energy window
            energy_range: record range of energy.
    */
    int i;
    double min_energy= *min_element(energy_distribution.begin(),energy_distribution.end());
    double max_energy= *max_element(energy_distribution.begin(),energy_distribution.end());
    int block_number= 20;  // use to specify number of energy blocks we want to consider.
    double energy_step = (max_energy- min_energy)/block_number;
    int state_number = energy_distribution.size();
    int index;
    state_density.assign(block_number,0);
    for(i=0;i<state_number;i++){
        index = (energy_distribution[i] - min_energy) / energy_step;
        if(index == block_number ) index= index -1;  // in case we meet max_energy element
        state_density[index]++;
    }
    double energy=min_energy;
    for(i=0;i<block_number+1;i++){
        energy_range.push_back(energy);
        energy= energy + energy_step;
    }
}

void detector:: output_state_density(vector<double> & dmat0,  vector<double> & dmat1){
    int i;
    // compute_detector state_density
    vector<int> state_density0;
    vector <int> state_density1;
    vector <double> energy_range0;
    vector<double> energy_range1;
    int block_number;

    vector<double> dmat_energy_level0;
    vector<double> dmat_energy_level1;
    if(my_id == 0) {
        ofstream state_density_output(path + "detector_state_density");
        for (i = 0; i < total_dmat_size[0]; i++) {
            dmat_energy_level0.push_back(dmat0[i]);  // dmat0 is diagonal part in all matrix.
        }
        for (i = 0; i < total_dmat_size[1]; i++) {
            dmat_energy_level1.push_back(dmat1[i]);
        }
        compute_detector_state_density(state_density0, dmat_energy_level0, energy_range0);
        compute_detector_state_density(state_density1, dmat_energy_level1, energy_range1);
        block_number = state_density0.size();
        state_density_output << "Detector 1 Range of energy block: " << endl;
        for (i = 0; i <= block_number; i++) {
            state_density_output << energy_range0[i] << " ";
        }
        state_density_output << endl;
        state_density_output << "Detector 1 density of state " << endl;
        for (i = 0; i < block_number; i++) {
            state_density_output << state_density0[i] << " ";
        }
        state_density_output << endl;
        state_density_output << " Detector 2 Range of energy block" << endl;
        for (i = 0; i <= block_number; i++) {
            state_density_output << energy_range1[i] << " ";
        }
        state_density_output << endl;
        state_density_output << "Detector 2 density of state " << endl;
        for (i = 0; i < block_number; i++) {
            state_density_output << state_density1[i] << " ";
        }
        state_density_output << endl;
        // output initial state's energy
        for (i = 0; i < 2; i++) {
            state_density_output << initial_state_energy[i] << "  " ;
        }
        state_density_output.close();
    }
}

void detector:: compute_local_density_of_state(ofstream & local_state_density_output,vector<double> & dmat0 , vector<double> & dmat1  ){
    // using eq.(2) in https://doi.org/10.1063/1.476070 to compute density of states: sum Lij
    // compute distribution of local_state_density
    // dmat0 contain all state's energy across process for detector0 . dmat1 contain all state's energy across process for detector 1
    int i,j, k ;
    double energy_difference;
    double ** density_of_states = new double * [electronic_state_num];
    double ** total_density_of_states = new double * [electronic_state_num];

    int ** dmatsize_displacement = new int * [electronic_state_num];
    for(i=0; i < electronic_state_num; i++){
        dmatsize_displacement[i] = new int [num_proc];
        dmatsize_displacement[i][0] = 0;
        for(j=1;j<num_proc;j++){
            dmatsize_displacement[i][j] = dmatsize_displacement[i][j-1] + dmatsize_each_process[i][j-1];
        }
    }
    int * begin_index = new int [2];

    int local_state_index = 0;

    int global_state_index_dirow = 0;
    int global_state_index_dicol = 0;

    double * coupling_strength_sum = new double [electronic_state_num];
    double * coupling_strength_sum_all_process = new double[electronic_state_num];
    double * coupling_strength_average_all_process = new double[electronic_state_num];
    vector<double> * dmat_ptr ;

    for(i=0; i < electronic_state_num; i++){
        begin_index[i] = total_dmat_size[i]/num_proc * my_id ;
    }

    for(i=0; i < electronic_state_num; i++){
        density_of_states[i] = new double [dmatsize[i]];
        for (j=0; j<dmatsize[i] ; j++ ){
            density_of_states[i][j] = 0;
        }
    }

    if(my_id == 0 ){
        for(i=0; i < electronic_state_num; i++){
            total_density_of_states[i] = new double [total_dmat_size[i]];
        }
    }

    for(i=0; i < electronic_state_num; i++){

        if(i==0){
            dmat_ptr = & dmat0;
        }
        else{
            dmat_ptr = & dmat1;
        }

        for( j = dmatsize[i]; j<dmatnum[i]; j++ ){
            local_state_index = dirow[i][j] - begin_index[i];
            global_state_index_dirow = dirow[i][j];
            global_state_index_dicol = dicol[i][j];

            energy_difference = abs ( (*dmat_ptr)[global_state_index_dirow] - (*dmat_ptr)[global_state_index_dicol] );
            density_of_states[i][local_state_index] = density_of_states[i][local_state_index] + 1 / (1 + pow(energy_difference / dmat[i][j],2) );
        }
    }
    // now we collect density_of_states in each process to total_density_of_states
    for(i=0; i < electronic_state_num; i++){
        MPI_Gatherv(&density_of_states[i][0],dmatsize[i],MPI_DOUBLE,&total_density_of_states[i][0],
                    dmatsize_each_process[i],dmatsize_displacement[i],MPI_DOUBLE,0,MPI_COMM_WORLD);
    }

    // compute coupling strength
    for(i=0; i < electronic_state_num; i++){
        coupling_strength_sum[i] = 0;
        for(j= dmatsize[i] ; j<dmatnum[i] ; j++){
            coupling_strength_sum [i] = coupling_strength_sum [i] + abs(dmat[i][j]);
        }
    }

    MPI_Allreduce(&coupling_strength_sum[0], &coupling_strength_sum_all_process[0], electronic_state_num, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    for(i=0; i < electronic_state_num; i++){
        coupling_strength_average_all_process[i] = coupling_strength_sum_all_process[i] / (total_dmat_num[i] - total_dmat_size[i]);
    }

    if(my_id == 0){
        // output result
        // format : number of states.  quantum number of state | 0 2 1 0 0 0 >. local density of state for that state.
        for(i=0; i < electronic_state_num; i++){
            local_state_density_output << coupling_strength_average_all_process[i] << endl;
        }

        for(i=0; i < electronic_state_num; i++){
            local_state_density_output << total_dmat_size[i] << endl;
            for(j=0;j<total_dmat_size[i];j++){
                for(k=0;k<nmodes[i];k++){
                    local_state_density_output << dv_all[i][j][k] << " ";
                }
                local_state_density_output << endl;
                local_state_density_output << total_density_of_states[i][j] << endl;
            }
        }
    }

    delete [] begin_index;
    for(i=0; i < electronic_state_num; i++){
        delete [] density_of_states[i];
    }

    if(my_id == 0){
        for(i=0; i < electronic_state_num; i++){
            delete [] total_density_of_states[i];
        }
    }
    delete [] density_of_states;
    delete [] total_density_of_states;
    for(i=0; i < electronic_state_num; i++){
        delete [] dmatsize_displacement [i];
    }
    delete [] dmatsize_displacement;
    delete [] coupling_strength_sum;
    delete [] coupling_strength_sum_all_process;
    delete [] coupling_strength_average_all_process;
}