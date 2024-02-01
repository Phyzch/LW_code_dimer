//
// Created by phyzch on 7/13/20.
//
#include "../system.h"
#include"../util.h"
using namespace std;



void full_system:: etot_MPI(double * hx, double * hy){
    int i;
    int irow_index, icol_index;
    double local_test=0;
    double total_test=0;
    double local_e=0;
    double local_e_add = 0;
    norm = 0 ;
    total_norm=0;
    // update tosend_xd, tosend_yd to send data to other process
    update_x_y();

    for(i=0;i<matsize;i++){
        hx[i]=0;
        hy[i]=0;
    }

    for(i=0;i<matnum;i++){
        irow_index = local_irow[i];
        icol_index = local_icol[i];
        hx[irow_index] = hx[irow_index] + mat[i] * real_part_wave_func[icol_index];
        hy[irow_index] = hy[irow_index] + mat[i] * imag_part_wave_func[icol_index];
    }

    for(i=0;i<matsize;i++){
        local_e = local_e + real_part_wave_func[i] * hx[i] + imag_part_wave_func[i] * hy[i];
        norm = norm + real_part_wave_func[i] * real_part_wave_func[i] + imag_part_wave_func[i] * imag_part_wave_func[i];
        local_test= local_test + hx[i] * imag_part_wave_func[i] - hy[i] * real_part_wave_func[i];
    }
    total_energy=0;
    MPI_Allreduce(&local_e,&total_energy,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&norm,&total_norm,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&local_test,&total_test,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    if(abs(total_test)>0.01){
        if(my_id==0){
            log << "Warning: The Hermiticity  is not conserved in Etot " << endl;
            cout << "Warning: The Hermiticity  is not conserved in Etot " << endl;
            MPI_Abort(MPI_COMM_WORLD,-1);
        }
    }
    if(abs(total_norm -1)>0.01){
        cout<<"Wrong"<<endl;
        if(my_id == 0){
            log <<" Warning the norm of whole system is not conserved "<<endl;
            cout <<"Warning the norm of whole system is not conserved"<<endl;
            cout<<total_norm<<endl;
        }
    }
    total_energy= total_energy/ total_norm;
}

void full_system :: generate_label_for_electronic_survival_prob_calculation(vector<double> & electronic_state_label_array){
    // for states with electronic index = 0, we set electronic state label = 1.
    // for states with electronic index = 1, we set electronic state label = 0.
    // therefore we can compute electronic state survival probability.
    int i,j, k;
    double electronic_state_label;
    for(i=0;i<matsize;i++){
        if(exciton_state_index_list[i] == 0){
            electronic_state_label = 1;
        }
        else{
            electronic_state_label = 0;
        }
        electronic_state_label_array.push_back(electronic_state_label);
    }

}


