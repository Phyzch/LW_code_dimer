// Chenghao Zhang  cz38@illinois.edu
#include"system.h"
#include"util.h"
using namespace std;
namespace fs = std::experimental::filesystem;

// Information for MPI program
int my_id;
int num_proc;
bool self_anharmonicity_bool = true;
double self_anharmonicity_D = 30000;


void output_survival_prob( string file_path , vector<vector<vector<int>>> & state_quantum_number_list,
                           vector<double> & state_energy_list,
                           vector<vector<double>> & time_list_all_states, vector<vector<double>> & survival_prob_list_all_states,
                           vector<vector<double>> & electronic_survival_prob_list_all_states );

void output_Nloc_for_vibrational_dimer_states(string file_path , vector<vector<vector<int>>> & state_quantum_number_list,
                                              vector<double> & state_energy_list,
                                              vector<vector<vector<int>>> & coupling_state_index_list,
                                              vector<vector<vector<vector<int>>>> & coupling_state_qn_list,
                                              vector<vector<vector<double>>> & coupling_state_strength_list,
                                              vector<vector<vector<double>>> & coupling_state_energy_diff_list,
                                              vector<vector<double>> & effective_coupling_number_list);

void output_vibrational_energy_dimer_states(string file_path, vector<vector<vector<int>>> & state_quantum_number_list,
                                            vector<double> & state_energy_list,
                                            vector<vector<double>> & time_list_all_states,
                                            vector<vector<vector<double>>> & vibrational_energy_list_all_states);

void read_state_quantum_number_list(string file_path, vector<vector<vector<int>>> & state_quantum_number_list);

int main(int argc,char * argv []) {
    srand(time(0));
    string path = "//home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/BChl_dimer_model/5_mode/"
                  "batch_simulation_Bigwood_scaling/nonstatistical_states/try/";


    string s;
    string s1;
    int i,m,j;
    int state_num;
    int mode_num;

    // MPI Command
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,&num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);

    vector<vector<int>> state_quantum_number;
    vector<vector<vector<int>>> state_quantum_number_list;

    read_state_quantum_number_list(path, state_quantum_number_list);
    state_num = state_quantum_number_list.size();
    mode_num = state_quantum_number_list[0][0].size();

    vector<double> state_energy_list;
    vector<vector<double>> time_list_all_states;
    vector<vector<double>> survival_prob_list_all_states;
    vector<vector<double>> electronic_survival_prob_list_all_states;

    // for local density of state Nloc:
    vector<vector<vector<int>>>  coupling_state_index_list;
    vector<vector<vector<vector<int>>>>  coupling_state_qn_list;
    vector<vector<vector<double>>>  coupling_state_strength_list;
    vector<vector<vector<double>>>  coupling_state_energy_diff_list;
    vector<vector<double>>  effective_coupling_number_list;

    // for vibrational energy in each monomer
    vector<vector<vector<double>>> vibrational_energy_list_all_states;

    if(my_id == 0){
        cout << "total state num for simulation : " << state_num << endl;
    }
    for(i=0; i < state_num; i++){

        { // the parenthese here let destructor called after we use this instance.
            state_quantum_number = state_quantum_number_list[i];

            vector<double> time_list;
            vector<double> survival_prob_list;
            vector<double> electronic_survival_prob_list;
            double state_energy;

            vector<vector<double>> vibrational_energy_list;

            full_system photon_entangled_system(path , state_quantum_number);  // set parameter and construct Hamiltonian.

            // compute Nloc for each state.
            photon_entangled_system.compute_local_density_of_state(coupling_state_index_list, coupling_state_qn_list,
                                                                   coupling_state_strength_list, coupling_state_energy_diff_list, effective_coupling_number_list);


            photon_entangled_system.Quantum_evolution(state_energy, time_list, survival_prob_list, electronic_survival_prob_list,
                                                      vibrational_energy_list); // creat initial state (or read from file). Then complete simulation.


            state_energy_list.push_back(state_energy);
            time_list_all_states.push_back(time_list);
            survival_prob_list_all_states.push_back(survival_prob_list);
            electronic_survival_prob_list_all_states.push_back(electronic_survival_prob_list);
            vibrational_energy_list_all_states.push_back(vibrational_energy_list);
        }

        if(my_id == 0){
            cout << "finish state num " << i << "   quantum number: " ;
            for(m=0;m<2;m++){
                for(j=0;j<mode_num;j++){
                    cout << state_quantum_number[m][j] << " ";
                }
                cout << "       " ;
            }
            cout << endl;
        }

    }

    output_survival_prob(path, state_quantum_number_list, state_energy_list, time_list_all_states, survival_prob_list_all_states,
                         electronic_survival_prob_list_all_states);

    output_Nloc_for_vibrational_dimer_states( path , state_quantum_number_list,
            state_energy_list, coupling_state_index_list,
            coupling_state_qn_list, coupling_state_strength_list,
            coupling_state_energy_diff_list, effective_coupling_number_list);

    output_vibrational_energy_dimer_states(path, state_quantum_number_list, state_energy_list, time_list_all_states,
                                           vibrational_energy_list_all_states);

    MPI_Finalize();
}

void read_state_quantum_number_list(string file_path, vector<vector<vector<int>>> & state_quantum_number_list){
    // file_name : sampling_state_info.txt
    // structure : state_num , mode_num_each_monomer. Then mode num for dimer each line.
    int state_num;
    int mode_num_each_monomer;
    int i,j;
    int m;
    int quantum_number;

    // read state quantum number from the file
    ifstream state_input;
    state_input.open(file_path + "sampling_state_info.txt");
    if(! state_input.is_open() ){
        if(my_id == 0){
            cout<< "THE INPUT STATE FILE FAILS TO OPEN!" << endl;
            MPI_Abort(MPI_COMM_WORLD, -35);
        }

    }

    state_input >> state_num >> mode_num_each_monomer;
    for(i=0;i<state_num;i++){
        vector<vector<int>> dimer_qn;
        for(m=0;m<2;m++){
            vector<int> monomer_qn;
            for(j=0;j<mode_num_each_monomer;j++){
                state_input >> quantum_number;
                monomer_qn.push_back(quantum_number);
            }
            dimer_qn.push_back(monomer_qn);
        }
        state_quantum_number_list.push_back(dimer_qn);
    }

}

void output_survival_prob( string file_path , vector<vector<vector<int>>> & state_quantum_number_list,
                    vector<double> & state_energy_list,
                    vector<vector<double>> & time_list_all_states, vector<vector<double>> & survival_prob_list_all_states,
                    vector<vector<double>> & electronic_survival_prob_list_all_states ){
    int i, j, k, m;
    int electronic_state_num = 2;
    int nmodes = state_quantum_number_list[0][0].size();
    vector<double> time_list = time_list_all_states[0];
    int time_step = time_list.size();

    ofstream survival_prob_out;
    ofstream electronic_survival_prob_out;

    int state_num = state_energy_list.size();

    if(my_id == 0){
        survival_prob_out.open(file_path + "survival_prob.txt");
        survival_prob_out << state_num << endl;
        // output state energy
        for(i=0; i< state_num; i++){
            survival_prob_out << state_energy_list[i] << " ";
        }
        survival_prob_out << endl;

        // output mode quanta
        for(i=0;i < state_num; i++){
            for(m=0; m< electronic_state_num; m++){
                for(k=0; k<nmodes; k++){
                    survival_prob_out << state_quantum_number_list[i][m][k] << " ";
                }
                survival_prob_out << "       ";
            }
            survival_prob_out << endl;
        }
        // output survival probability
        for(i=0;i<time_step;i++){
            survival_prob_out << time_list[i] << endl;
            for(j=0;j<state_num;j++){
                survival_prob_out << survival_prob_list_all_states[j][i] << " ";
            }
            survival_prob_out << endl;
        }
    }

    if(my_id == 0){
        electronic_survival_prob_out.open(file_path + "electronic_survival_prob.txt");
        electronic_survival_prob_out <<  state_num << endl;
        // output state energy
        for(i=0; i< state_num; i++){
            electronic_survival_prob_out << state_energy_list[i] << " ";
        }
        electronic_survival_prob_out << endl;
        // output mode quanta
        for(i=0;i < state_num; i++){
            for(m=0; m< electronic_state_num; m++){
                for(k=0; k<nmodes; k++){
                    electronic_survival_prob_out << state_quantum_number_list[i][m][k] << " ";
                }
                electronic_survival_prob_out << "       ";
            }
            electronic_survival_prob_out << endl;
        }

        // output electronic survival prob
        for(i=0;i < time_step;i++){
            electronic_survival_prob_out << time_list[i] << endl;
            for(j=0;j<state_num;j++){
                electronic_survival_prob_out << electronic_survival_prob_list_all_states[j][i] << " ";
            }
            electronic_survival_prob_out << endl;
        }

    }

}

void output_Nloc_for_vibrational_dimer_states(string file_path , vector<vector<vector<int>>> & state_quantum_number_list,
                                              vector<double> & state_energy_list,
                                              vector<vector<vector<int>>> & coupling_state_index_list,
                                              vector<vector<vector<vector<int>>>> & coupling_state_qn_list,
                                              vector<vector<vector<double>>> & coupling_state_strength_list,
                                              vector<vector<vector<double>>> & coupling_state_energy_diff_list,
                                              vector<vector<double>> & effective_coupling_number_list){
    int i, j, k, m;
    int electronic_state_num = 2;
    int nmodes = state_quantum_number_list[0][0].size();
    int state_num = state_quantum_number_list.size();
    int coupled_state_num;

    double coupling_strength;
    double energy_difference;
    double effective_number_contribution;

    double effective_coupling_number_dimer;

    string coupling_info_file_name = "same_PES_coupling_info.txt";
    string Nloc_file_name = "same_PES_Nloc.txt";
    coupling_info_file_name = file_path + coupling_info_file_name;
    Nloc_file_name = file_path + Nloc_file_name;

    ofstream coupling_info_output;
    ofstream Nloc_output;
    if(my_id == 0){
        coupling_info_output.open(coupling_info_file_name);
        coupling_info_output << state_num << "  " << nmodes << endl;
        coupling_info_output << "quantum_number" << "  " << "coupling_strength" << "  " << "energy_difference" << " Nloc_contribution" << endl;
        coupling_info_output << endl;

        for(i=0;i<state_num;i++){
            for(m=0;m<electronic_state_num;m++){
                for(j=0;j<nmodes;j++){
                    coupling_info_output << state_quantum_number_list[i][m][j] << " ";
                }
                coupling_info_output << "        ";
            }
            coupling_info_output << state_energy_list[i] << endl;
            for(m=0;m<electronic_state_num;m++){
                for(j=0;j<nmodes;j++){
                    coupling_info_output << state_quantum_number_list[i][m][j] << " ";
                }
                coupling_info_output << endl;
                coupling_info_output << effective_coupling_number_list[i][m] << endl;

                coupled_state_num = coupling_state_qn_list[i][m].size();
                for(k=0;k<coupled_state_num;k++){
                    for(j=0;j<nmodes;j++){
                        coupling_info_output << coupling_state_qn_list[i][m][k][j] << " ";
                    }
                    coupling_strength = coupling_state_strength_list[i][m][k];
                    energy_difference = coupling_state_energy_diff_list[i][m][k];
                    effective_number_contribution = 1 / (1 + pow(energy_difference / coupling_strength , 2));

                    coupling_info_output << "  " << coupling_strength << "  " << energy_difference << "  " <<  effective_number_contribution << endl;
                }
                coupling_info_output << endl;
            }
            coupling_info_output << endl;
        }
        coupling_info_output.close();
    }

    if(my_id == 0){
        Nloc_output.open(Nloc_file_name);
        Nloc_output << state_num << "  " << nmodes << endl;
        Nloc_output << "quantum_number" << "  " << "coupling_strength" << "  " << "energy_difference" << " Nloc_contribution" << endl;
        Nloc_output << endl;

        for(i=0;i<state_num;i++){
            for(m=0;m<electronic_state_num;m++){
                for(j=0;j<nmodes;j++){
                    Nloc_output << state_quantum_number_list[i][m][j] << " ";
                }
                Nloc_output << "        ";
            }
            Nloc_output << state_energy_list[i] << endl;

            effective_coupling_number_dimer = effective_coupling_number_list[i][0] + effective_coupling_number_list[i][1];
            Nloc_output << effective_coupling_number_dimer << endl;

            for(m=0;m<electronic_state_num;m++){
                for(j=0;j<nmodes;j++){
                    Nloc_output << state_quantum_number_list[i][m][j] << " ";
                }
                Nloc_output << endl;
                Nloc_output<< effective_coupling_number_list[i][m];
                Nloc_output << endl;
            }
            Nloc_output<< endl;
        }
        Nloc_output.close();
    }

}

void output_vibrational_energy_dimer_states(string file_path, vector<vector<vector<int>>> & state_quantum_number_list,
                                            vector<double> & state_energy_list,
                                            vector<vector<double>> & time_list_all_states,
                                            vector<vector<vector<double>>> & vibrational_energy_list_all_states){
    int i, j, k, m;
    int electronic_state_num = 2;
    int nmodes = state_quantum_number_list[0][0].size();
    vector<double> time_list = time_list_all_states[0];
    int time_step = time_list.size();

    ofstream vibrational_energy_output;
    int state_num = state_energy_list.size();

    string file_name = "vibrational_energy.txt";

    if(my_id == 0){
         vibrational_energy_output.open(file_path + file_name);
         vibrational_energy_output << state_num << endl;
        // output state energy
        for(i=0; i< state_num; i++){
            vibrational_energy_output << state_energy_list[i] << " ";
        }
        vibrational_energy_output << endl;
        // output mode quanta
        for(i=0;i < state_num; i++){
            for(m=0; m< electronic_state_num; m++){
                for(k=0; k<nmodes; k++){
                    vibrational_energy_output << state_quantum_number_list[i][m][k] << " ";
                }
                vibrational_energy_output << "       ";
            }
            vibrational_energy_output << endl;
        }
        // output vibrational energy
        for(i=0;i<time_step;i++){
            vibrational_energy_output << time_list[i] << endl;
            for(m=0;m<electronic_state_num;m++){
                for(j=0;j<state_num;j++){
                    vibrational_energy_output << vibrational_energy_list_all_states[j][m][i] << " ";
                }
                vibrational_energy_output << endl;
            }
        }
        vibrational_energy_output.close();
    }

}