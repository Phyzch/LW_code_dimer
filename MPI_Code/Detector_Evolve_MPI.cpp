//
// Created by phyzch on 6/27/20.
//
#include "../system.h"
#include "../util.h"
using namespace std;

int  detector::construct_receive_buffer_index(int * remoteVecCount_element, int * remoteVecPtr_element, int * remoteVecIndex_element, int detector_index){
    // input: remoteVecCount: total number of element need to receive from each process.
    //        remoteVecPtr: displacement in remoteVecIndex for element in each process.
    //        remoteVecIndex: index for remote vector we need to receive. (they may allocate in different remote process.)
    // return: length of remoteVecIndex.

    int i,j;
    // range for element in process is [local_begin, local_end)
    int total_remoteVecCount=0;
    int vsize = total_dmat_size[detector_index] / num_proc;
    int local_begin= total_dmat_size[detector_index] / num_proc * my_id;
    int local_end;
    int remote_pc_id;
    if(my_id!=num_proc-1) {
        local_end = total_dmat_size[detector_index] / num_proc * (my_id + 1);
    }
    else{
        local_end = total_dmat_size[detector_index];
    }
    // ---------------------------------------------------------------
    vector <int> col_index_copy = dicol[detector_index];
    sort(col_index_copy.begin(),col_index_copy.end()); // sort vector.
    int col_array_size = col_index_copy.size();
    int prev_col=-1;
    j=0;
    for(i=0;i<col_array_size;i++){
        if( (col_index_copy[i]>prev_col)   and ( (col_index_copy[i]<local_begin) or (col_index_copy[i] >= local_end) )  ){
            // this matrix element is not in process.
            if (col_index_copy[i] >= vsize * (num_proc-1) ){
                remote_pc_id = num_proc-1;
            }
            else{
                remote_pc_id = col_index_copy[i] / vsize;
            }
            remoteVecCount_element[remote_pc_id] ++;
            remoteVecIndex_element [j] = col_index_copy[i];  // vector index need to receive. (global index , ordered)
            j++;
        }
        prev_col= col_index_copy[i];
    }
    remoteVecPtr_element[0]=0;   // displacement for remote vector from each process in remoteVecIndex.
    for(i=1;i<num_proc;i++){
        remoteVecPtr_element[i] = remoteVecPtr_element[i-1] + remoteVecCount_element[i-1];
    }
    for(i=0;i<num_proc;i++){
        total_remoteVecCount = total_remoteVecCount + remoteVecCount_element[i];
    }
    return total_remoteVecCount;
}

int construct_send_buffer_index(int * remoteVecCount_element, int * remoteVecPtr_element, int * remoteVecIndex_element,
                                int * tosendVecCount_element, int * tosendVecPtr_element, int* & tosendVecIndex_ptr){
    //  tosend_Vec_count record number of element to send to each process.
    // tosend_Vec_Index record the global index of vector the process have to send
    // tosend_Vec_Ptr record the offset of vector to send to each other process.
    // return to_send_buffer_len: lenfth of tosendVecIndex
    int i;
    int to_send_buffer_len;

    MPI_Alltoall(&remoteVecCount_element[0],1,MPI_INT,&(tosendVecCount_element[0]),1,MPI_INT,MPI_COMM_WORLD);

    // compute displacement for each process's data.
    tosendVecPtr_element[0]=0;
    for(i=1;i<num_proc;i++){
        tosendVecPtr_element[i]= tosendVecPtr_element[i-1] + tosendVecCount_element[i-1];
    }
    // compute total length of buffer to send
    to_send_buffer_len=0;
    for(i=0;i<num_proc;i++){
        to_send_buffer_len= to_send_buffer_len + tosendVecCount_element[i];
    }
    // Index (in global) of element to send. use MPI_Alltoallv to receive the index to send.
    tosendVecIndex_ptr = new int [to_send_buffer_len];
    MPI_Alltoallv(&remoteVecIndex_element[0],remoteVecCount_element,remoteVecPtr_element,MPI_INT,
                  & tosendVecIndex_ptr[0],tosendVecCount_element,tosendVecPtr_element,MPI_INT,MPI_COMM_WORLD);

    return to_send_buffer_len;
}

int compar(const void * a, const void * b){
    return *(int *) a - * (int *)b;
}

// this function is called every time we do evolution
void detector::prepare_evolution(){
    // compute buffer to receive and send for each process.
    // resize xd,yd to provide extra space for recv_buffer.
    // allocate space for send_xd , send_yd buffer.
    // Index for remoteVecIndex, tosendVecIndex are computed here.
    int m,i;
    int vsize;

    // Index for vector to send and receive.
    // remoteVecCount: total number to receive. remoteVecPtr: displacement in remoteVecIndex for each process. remoteVecIndex: index in other process to receive.
    // tosendVecCount: total number to send to other process. tosendVecPtr: displacement in tosendVecIndex in each process.  tosendVecIndex: Index of element in itself to send. (it's global ,need to be converted to local index)
    remoteVecCount= new int * [stlnum];
    remoteVecPtr= new int * [stlnum];
    remoteVecIndex= new int * [stlnum];
    to_recv_buffer_len = new int  [stlnum];

    //------------------Allocate space for vector to receive ---------------------
    for (m=0;m<stlnum;m++){
        remoteVecCount[m] = new int [num_proc];
        remoteVecPtr[m] = new int [num_proc];
        remoteVecIndex[m] = new int [dmatnum[m]];
        for(i=0;i<num_proc;i++){
            remoteVecCount[m][i] = 0;
        }
    }
    tosendVecCount= new int *[stlnum];
    tosendVecPtr =  new int * [stlnum];
    tosendVecIndex = new int * [stlnum];
    to_send_buffer_len= new int [stlnum];
    for(m=0;m<stlnum;m++){
        tosendVecCount[m] = new int [num_proc];
        tosendVecPtr[m] = new int [num_proc];
    }

    int * search_Ind; // local variable, used for compute local_dicol;
    int col_index_to_search;
    // local column index used when we do H *x and H*y
    local_dirow= new vector<int> [2];
    local_dicol = new vector<int> [2]; // column index for computation in local matrix.
    // buffer to send and receive buffer to/from other process.
    recv_xd= new double * [stlnum];
    recv_yd= new double * [stlnum];
    send_xd= new double * [stlnum];
    send_yd = new double *[stlnum];
    for(m=0;m<stlnum;m++){
        vsize= total_dmat_size[m]/num_proc;
        to_recv_buffer_len[m] = construct_receive_buffer_index(remoteVecCount[m],remoteVecPtr[m],
                remoteVecIndex[m],m);  // construct buffer to receive.
        to_send_buffer_len[m]= construct_send_buffer_index(remoteVecCount[m],remoteVecPtr[m],remoteVecIndex[m],
                                                        tosendVecCount[m],tosendVecPtr[m], tosendVecIndex[m]);
        xd[m].resize(dmatsize[m] + to_recv_buffer_len[m]);
        yd[m].resize(dmatsize[m] + to_recv_buffer_len[m]);
        recv_xd[m] = new double [to_recv_buffer_len[m]];
        recv_yd[m]= new double [to_recv_buffer_len[m]];
        send_xd[m]= new double [to_send_buffer_len[m]];
        send_yd[m] = new double [to_send_buffer_len[m]];
        // construct local_dirow, local_dicol
        local_dicol[m].reserve(dmatnum[m]);
        local_dirow[m].reserve(dmatnum[m]);
        for(i=0;i<dmatnum[m];i++){
            local_dirow[m].push_back(dirow[m][i] - my_id * vsize);  // set local index for row index
            col_index_to_search= dicol[m][i];
            search_Ind=(int *) bsearch(&col_index_to_search,remoteVecIndex[m],to_recv_buffer_len[m],sizeof(int),compar);
            if(search_Ind!=NULL){
                // this column index is not in local matrix, and we should get it from other process (remoteVec)
                local_dicol[m].push_back(dmatsize[m] + (search_Ind-remoteVecIndex[m]) );
            }
            else{ // this column index is in local matrix.
                local_dicol[m].push_back (dicol[m][i] - my_id * vsize );
            }
        }
    }

}
void detector::update_dx_dy(int detector_index){
    int i,m;
    int vsize;
    m= detector_index;
    // collect data for send_buffer.
    vsize = total_dmat_size[m]/num_proc;
    for (i = 0; i < to_send_buffer_len[m]; i++) {
        send_xd[m][i] = xd[m][tosendVecIndex[m][i] - my_id * vsize];
        send_yd[m][i] = yd[m][tosendVecIndex[m][i] - my_id * vsize];
    }

    MPI_Alltoallv(&send_xd[m][0],tosendVecCount[m],tosendVecPtr[m],MPI_DOUBLE,
                  &recv_xd[m][0],remoteVecCount[m],remoteVecPtr[m],MPI_DOUBLE,MPI_COMM_WORLD);
    MPI_Alltoallv(&send_yd[m][0],tosendVecCount[m],tosendVecPtr[m],MPI_DOUBLE,
                  &recv_yd[m][0],remoteVecCount[m],remoteVecPtr[m],MPI_DOUBLE,MPI_COMM_WORLD);

    // copy received array to xd, yd:  This part may be optimized for example directly send result in xd, yd.
    for(i=0;i<to_recv_buffer_len[m];i++){
        xd[m][i+ dmatsize[m]]= recv_xd[m][i];
        yd[m][i+ dmatsize[m]]= recv_yd[m][i];
    }

}
// call it every time we do SUR algorithm. do it for detector specified by detector_index.
void detector::SUR_onestep_MPI(int detector_index, double cf){
    int m,i;
    int irow,icol;
    m = detector_index;
    update_dx_dy(detector_index);
    for(i=0;i<dmatnum[m];i++){
        // make sure when we compute off-diagonal matrix, we record both symmetric and asymmetric part
        irow = local_dirow[m][i];
        icol = local_dicol[m][i]; // compute to point to colindex in
        xd[m][irow] = xd[m][irow] + dmat[m][i] * yd[m][icol] * cf;
    }

    update_dx_dy(detector_index);

    for(i=0;i<dmatnum[m];i++){
        irow= local_dirow[m][i];
        icol= local_dicol[m][i];
        yd[m][irow] = yd[m][irow] - dmat[m][i] * xd[m][icol] * cf;
    }
}

void full_system::pre_coupling_evolution_MPI(int initial_state_choice){
    // we do not write output function now, try to make it as simple as possible.
    int irow_index, icol_index;
    int m,i,j,k;
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    int steps;
    double detector_tprint = 0.01;
    int output_step= int(detector_tprint/delt); //Output every output_step.
    double * start_time = new double [s.tlnum];
    double * final_time = new double [s.tlnum];
    for(i=0;i<s.tlnum;i++){
        start_time[i]=0;
    }
    ofstream Detector_precoup_mode_quanta;
    ofstream Detector_precoup_output;
    ofstream Detector_energy;
    // -------------used for find special state (bright state or initial state ) to print out survival probability
    int ** special_state;
    int * special_state_pc_id, *special_state_index;
    double * special_state_x, * special_state_y;
    bool exist;
    int position;
    // -----------Open detector_precoup_mode_quanta ofstream -----------------------------
    if(my_id==0){
        if(initial_state_choice==1) {
            Detector_precoup_mode_quanta.open(path + "bright_state_detector_precoup_mode_quanta.txt");
            Detector_precoup_output.open(path + "bright_state_detector_precoupling.txt");
            Detector_energy.open(path+"bright_state_detector_energy.txt");
        }
        else{
            Detector_precoup_mode_quanta.open(path + "lower_bright_state_detector_precoup_mode_quanta.txt");
            Detector_precoup_output.open(path + "lower_bright_state_detector_precoupling.txt");
            Detector_energy.open(path+"lower_bright_state_detector_energy.txt");
        }
    }

    //---------- Allocate space for mode quanta -------------------------------------------
    double ** total_mode_quanta;
    total_mode_quanta= new double * [s.tlnum];
    for(m=0;m<s.tlnum;m++){
        total_mode_quanta[m]= new double [d.nmodes[m]];
    }

    double ** mode_quanta= new double * [s.tlnum];
    for (m=0;m<s.tlnum;m++){
        mode_quanta[m]= new double [d.nmodes[m]];
    }
    //--------------Prepare to output detector bright/initial state  -----------------
    if(initial_state_choice==1){
        special_state= d.bright_state;
    }
    else{
        special_state=d.initial_detector_state;
    }
    special_state_pc_id = new int [s.tlnum];  // record process id that record special state wave function.
    special_state_index = new int [s.tlnum]; // index for special state in process: special_state_pc_id.
    special_state_x = new double [s.tlnum];
    special_state_y= new double [s.tlnum];
    bool * exist_bool_for_pc = new bool [num_proc];
    for(m=0;m<s.tlnum;m++){
        vector <int> vec_special_state;
        for(i=0;i<d.nmodes[m]; i++){
            vec_special_state.push_back(special_state[m][i]);
        }
        position=find_position_for_insert_binary(d.dv[m],vec_special_state,exist);
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
                log<<" Can not find initial state or brigth state in all vmode. Must have bug here."<<endl;
                MPI_Abort(MPI_COMM_WORLD,-7);
            }
        }
        // tell every process the location of special state
        MPI_Bcast(&position,1, MPI_INT,special_state_pc_id[m],MPI_COMM_WORLD);
        special_state_index [m] = position;
    }
    // -----------------------------------------------------------------------------------------------
    // prepare sendbuffer and recv_buffer and corresponding index.
    d.prepare_evolution();
    vector<complex<double>> * H_phi = new vector<complex<double>> [s.tlnum];
    double de;
    double de_all;
    for(m=0;m<s.tlnum;m++) {
        H_phi[m].resize(d.dmatsize[m]);
        for(i=0;i<d.dmatsize[m];i++){
            H_phi[m][i] = 0;
        }
    }
    for(m=0;m<s.tlnum;m++){
        final_time[m]=0;
        if(d.proptime[m]>0){
            if(my_id==0){
                log <<" Pre-propogation for detector  "<<m<<endl;
            }
            t=0;
            steps= d.proptime[m]/delt + 1;
            if(my_id==0){
                Detector_precoup_mode_quanta << "Detector mode quanta for Detector " << m << endl;
                Detector_precoup_mode_quanta << "total time: " << (d.proptime[m]-start_time[m]) << " " << delt * output_step << endl;
                Detector_precoup_output << "Wavefunction for Detector  " << m << endl;
                Detector_precoup_output << "total time: " << (d.proptime[m]-start_time[m]) << " " << detector_tprint
                                        << endl;
            }
            // Do simulation in loop
            for(k=0;k<steps;k++){
                //-------------------- output result ----------------------------
                if(k % output_step ==0) {
                    // ---------------------------------output detector mode quanta -------------------------------------------------
                    for (j = 0; j < d.nmodes[m]; j++) {
                        mode_quanta[m][j] = 0;
                    }
                    for (i = 0; i < d.dmatsize[m]; i++) {
                        for (j = 0; j < d.nmodes[m]; j++) {
                            mode_quanta[m][j] =
                                    mode_quanta[m][j] + (pow(d.xd[m][i], 2) + pow(d.yd[m][i], 2)) * d.dv[m][i][j];
                        }
                    }
                    MPI_Reduce(&mode_quanta[m][0], &total_mode_quanta[m][0], d.nmodes[m], MPI_DOUBLE, MPI_SUM, 0,
                               MPI_COMM_WORLD);
                    if (my_id == 0) {
                        Detector_precoup_mode_quanta << "Time:   " << t << endl;
                        for (j = 0; j < d.nmodes[m]; j++) {
                            Detector_precoup_mode_quanta << total_mode_quanta[m][j] << " ";
                        }
                        Detector_precoup_mode_quanta << endl;
                    }

                    //---------------------- output detector special state (excited bright state or lower bright state)-
                    if (my_id == special_state_pc_id[m]) {
                            special_state_x[m] = d.xd[m][special_state_index[m]];
                            special_state_y[m] = d.yd[m][special_state_index[m]];
                    }
                    MPI_Bcast(&special_state_x[m],1,MPI_DOUBLE,special_state_pc_id[m],MPI_COMM_WORLD);
                    MPI_Bcast(&special_state_y[m], 1, MPI_DOUBLE, special_state_pc_id[m], MPI_COMM_WORLD);
                    if (my_id == 0) {
                        // output result to Detector_wavefunction.txt
                        Detector_precoup_output << "Time:   " << t << endl;
                        Detector_precoup_output << "real part: ";
                        Detector_precoup_output << special_state_x[m] << endl;
                        Detector_precoup_output << "imaginary part: ";
                        Detector_precoup_output << special_state_y[m] << endl;
                    }
                    //--------------------------------------------------------------------------------------------------
                    // compute detector energy
                    for(i=0;i<d.dmatsize[m];i++) {
                        H_phi[m][i]=0;
                    }
                    for(i=0;i<d.dmatnum[m];i++){
                        irow_index = d.local_dirow[m][i];
                        icol_index = d.local_dicol[m][i]; // compute to point to colindex in
                        H_phi[m][irow_index] = H_phi[m][irow_index] + d.dmat[m][i] * complex(d.xd[m][icol_index],d.yd[m][icol_index]);
                    }
                    de=0;
                    for(i=0;i<d.dmatsize[m];i++){
                        de= de+ real(H_phi[m][i] * complex(d.xd[m][i],-d.yd[m][i]));
                    }
                    MPI_Allreduce(&de,&de_all,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
                    if(my_id == 0){
                        Detector_energy << de_all <<" ";
                    }


                }
                t= t+ delt;
                d.SUR_onestep_MPI(m,cf);
            }
            final_time[m]= t;
            Detector_energy<<endl;
        }


    }

//    d.save_detector_state_MPI(path,final_time,log,initial_state_choice);


    if(my_id==0){
        cout<<"Detector_pre_coupling simulation finished"<<endl;
        Detector_precoup_mode_quanta.close();
        Detector_precoup_output.close();
        Detector_energy.close();
    }
    // -------------- free remote_Vec_Count, remote_Vec_Index -------------------------
    for(i=0;i<s.tlnum;i++){
        delete [] d.remoteVecCount[i];
        delete [] d.remoteVecPtr[i];
        delete []  d.remoteVecIndex[i];
        delete [] d.tosendVecCount[i];
        delete [] d.tosendVecPtr[i];
        delete [] d.tosendVecIndex[i];
        delete [] d.send_xd[i];
        delete [] d.send_yd[i];
        delete [] d.recv_xd[i];
        delete [] d.recv_yd[i];
        delete [] total_mode_quanta[i];
        delete [] mode_quanta[i];
    }
    delete [] total_mode_quanta;
    delete [] mode_quanta;

    delete [] special_state_pc_id;
    delete [] special_state_index;
    delete [] special_state_x;
    delete [] special_state_y;

    delete [] exist_bool_for_pc;

    delete [] d.to_recv_buffer_len;
    delete [] d.remoteVecCount;
    delete [] d.remoteVecPtr;
    delete [] d.remoteVecIndex;
    delete [] d.to_send_buffer_len;
    delete [] d.tosendVecCount;
    delete [] d.tosendVecPtr;
    delete [] d.tosendVecIndex;
    delete [] d.send_xd;
    delete [] d.send_yd;
    delete [] d.recv_xd;
    delete [] d.recv_yd;
    delete [] d.local_dirow;
    delete [] d.local_dicol;

    if(initial_state_choice ==1){
        // if we simulate bright state, we will delete x[] ,y [] in case we may want to simulate detector lower bright state later.
        delete [] d.xd;
        delete [] d.yd;
    }
};

