//
// Created by phyzch on 2/10/23.
// compute reduced density matrix on 2 vibrational dof (exciton vibrational density)
//
#include "util.h"
#include "system.h"

void full_system::compute_exciton_vibrational_density( vector<vector<double>> & exciton_vib_density , int electronic_state ){
    // compute reduced density matrix in each electronic state.
    // the vibrational mode we choose is in monomer 1, index as specified in the code. (supposed to be 0 , 1)
    // electronic state indicate the electronic surface we choose to compute EVD for. can only be 0 or 1.

    int i,j,k;
    int mapping_mode0_qn, mapping_mode1_qn;
    int state_d1_index; // index for vibrational state in each monomer
    // initialize exciton vibrational density
    exciton_vib_density.clear();
    // used to compute exciton vibrational density for states in each process.
    vector<vector<double>> exciton_vib_density_each_pc;
    for(i=0; i<= d.nmax[0][mapping_mode0]; i++){
        vector<double> evd_row;
        for(j=0; j<= d.nmax[0][mapping_mode1]; j++){
            evd_row.push_back(0);
        }
        exciton_vib_density.push_back(evd_row);
        exciton_vib_density_each_pc.push_back(evd_row);
    }

    for(i=0; i<matsize;i++){
        if(sstate[i] != electronic_state) continue;
        state_d1_index = dstate[0][i];
        mapping_mode0_qn = d.dv_all[0][state_d1_index][mapping_mode0];
        mapping_mode1_qn = d.dv_all[0][state_d1_index][mapping_mode1];
        exciton_vib_density_each_pc[mapping_mode0_qn][mapping_mode1_qn] =
                exciton_vib_density_each_pc[mapping_mode0_qn][mapping_mode1_qn] + pow(x[i],2) + pow(y[i],2);
    }

    for(i=0;i<=d.nmax[0][mapping_mode0]; i++){
        MPI_Allreduce(&exciton_vib_density_each_pc[i][0], &exciton_vib_density[i][0], d.nmax[0][mapping_mode1] + 1, MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD);
    }

}