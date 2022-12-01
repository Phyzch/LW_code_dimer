//
// Created by phyzch on 11/30/22.
//

#include "../util.h"
#include "../system.h"
#include "../quotient_state.h"

int full_system::search_full_sys_matrix_given_sd_matrix(int s_state, int d1_state, int d2_state){
    // s_state : index for electronic state
    // d_state: index for detector1
    // d_state: index for detector2

    int x_index1; // index in full matrix.

    vector <int> d1_mode;
    vector <int> d2_mode;

    d1_mode = d.dv_all[0][d1_state];
    d2_mode = d.dv_all[1][d2_state];

    bool exist;
    int d_list_position;
    int xlist_position;

    // search from d1_list
    d_list_position = find_location_binarysearch_quotient_state(d1list, s_state, d2_mode, exist);
    if( not exist ){
        // this state does not exist
        return -1;
    }

    xlist_position = binary_search_dxindex(d1list[d_list_position].dxindex, d1_state, exist);
    if (not exist){
        return -1;
    }

    x_index1 = d1list[d_list_position].xindex[xlist_position];

    return x_index1;
}

double compute_franck_condon_factor(double alpha, int m, int n){
    // compute franck_condon factor <m|alpha;n>. Here alpha is displacement operator.
    int nm_min = min(n,m);
    int l;
    double franck_condon_factor = 0;
    // prefacto = e^{-alpha^2/2} sqrt(n! m!) * alpha^{n+m}
    // std::tgamma(n+1) = n!
    double prefactor = exp(- pow(alpha,2)/2) * pow(alpha, n + m) * sqrt(std::tgamma(n+1)  * std::tgamma(m+1) );

    if (alpha!=0){
        for (l=0;l<=nm_min;l++){
            // sum: 1/(l! * (n-l)! * (m-l)!) * (-1)^{n-l} * alpha^{-2l}
            franck_condon_factor = franck_condon_factor + 1/( std::tgamma(l+1) * std::tgamma(n-l+1) * std::tgamma(m-l+1) ) * pow(-1, n-l) * pow(alpha, -2*l);
        }

        franck_condon_factor = franck_condon_factor * prefactor;
    }

    if (alpha ==0){
        if (m==n){
            franck_condon_factor = 1;
        }
        else{
            franck_condon_factor = 0;
        }
    }

    return franck_condon_factor;
}


void detector::compute_franck_condon_factor_table() {
    // compute list of franck condon factor for nonadiabatic coupling.
    // different dof have different displacement (alpha).
    // results store in franck_condon_factor_table

    // max quantum number
    int max_qn = 0;
    int m;
    int i, j, k;
    double franck_condon_factor;
    double alpha;
    for(i=0;i<nmodes[0];i++){
        if (max_qn < nmax[0][i]) {
            max_qn = nmax[0][i];
        }
    }
    // include states from 0 to max_qn;
    max_qn = max_qn + 1;

    franck_condon_factor_table = new double *** [electronic_state_num];
    for (m = 0; m < electronic_state_num; m++){
        franck_condon_factor_table[m] = new double **[nmodes[m]];
        for(j=0;j<nmodes[m];j++){
            franck_condon_factor_table[m][j] = new double * [max_qn];
            for(k=0; k<nmodes[m]; k++){
                franck_condon_factor_table[m][j][k] = new double [max_qn];
            }
        }
    }

    for (m=0; m < electronic_state_num ; m++ ){
        for(i=0;i<nmodes[m];i++){
            alpha = electron_phonon_coupling[m][i];
            for(j=0;j<max_qn;j++){
                for(k=0;k<max_qn;k++){
                    franck_condon_factor = compute_franck_condon_factor(alpha, j , k);
                    franck_condon_factor_table[m][i][j][k] = franck_condon_factor;
                }
            }
        }
    }


    // output franck condon factor table <m| alpha;n>
    if (my_id == 0){
        ofstream franck_condon_output;
        franck_condon_output.open(path + "franck_condon.txt");
        for(m=0;m<electronic_state_num;m++){
            for(i=0;i<nmodes[m];i++){
                franck_condon_output << "alpha: " << electron_phonon_coupling[m][i] << endl;
                for(j=0;j<max_qn;j++){
                    for(k=0;k<max_qn;k++){
                        franck_condon_output << franck_condon_factor_table[m][i][j][k] << " ";
                    }
                    franck_condon_output << endl;
                }
                franck_condon_output << endl;
            }

            franck_condon_output << endl;
        }

    }

}

void detector::find_franck_condon_factor_for_monomer_states(){
    int i, j, k, l, m;
    double franck_condon_factor_prod;
    double franck_condon_factor;

    int state1_mode_qn, state2_mode_qn;

    // brute force method to find all d states
    for(m=0; m < electronic_state_num ; m++){
        vector<vector<int>> nonadiabatic_coupled_state_each_monomer;
        vector<vector<double>> nonadiabatic_coupled_state_each_monomer_fc;

        for(i=0; i< total_dmat_size[m] ; i++ ){
            vector<int> nonadiabatic_coupled_state_single_state;
            vector<double> nonadiabatic_coupled_state_single_state_fc;

            for( j=0; j<total_dmat_size[m]; j++ ){
                franck_condon_factor_prod = 1;

                for(k=0; k<nmodes[m]; k++){
                    state1_mode_qn = dv_all[m][i][k];
                    state2_mode_qn = dv_all[m][j][k];

                    franck_condon_factor = franck_condon_factor_table[m][k][state1_mode_qn][state2_mode_qn];
                    franck_condon_factor_prod = franck_condon_factor_prod * franck_condon_factor;
                }

                if( franck_condon_factor_prod > Franck_condon_factor_cutoff ){
                    nonadiabatic_coupled_state_single_state.push_back(j);
                    nonadiabatic_coupled_state_single_state_fc.push_back(franck_condon_factor_prod);
                }

            }
            nonadiabatic_coupled_state_each_monomer.push_back(nonadiabatic_coupled_state_single_state);
            nonadiabatic_coupled_state_each_monomer_fc.push_back(nonadiabatic_coupled_state_single_state_fc);
        }
        nonadiabatic_coupled_d_state.push_back(nonadiabatic_coupled_state_each_monomer);
        nonadiabatic_coupled_d_state_franck_condon.push_back(nonadiabatic_coupled_state_each_monomer_fc);

    }

}