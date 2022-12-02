// Chenghao Zhang  cz38@illinois.edu
#include"system.h"
#include"util.h"
using namespace std;
namespace fs = std::experimental::filesystem;

// Information for MPI program
int my_id;
int num_proc;

void output_survival_prob( string file_path , vector<vector<vector<int>>> & state_quantum_number_list,
                           vector<double> & state_energy_list,
                           vector<vector<double>> & time_list_all_states, vector<vector<double>> & survival_prob_list_all_states,
                           vector<vector<double>> & electronic_survival_prob_list_all_states );

void read_state_quantum_number_list(string file_path, vector<vector<vector<int>>> & state_quantum_number_list);

int main(int argc,char * argv []) {
    srand(time(0));
    string path = "/home/phyzch/Presentation/LW_electronic_model/2022 result/spin_boson_LW/BChl_dimer_model/5_mode/batch_simulation/try/";


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

            full_system photon_entangled_system(path , state_quantum_number);  // set parameter and construct Hamiltonian.
            photon_entangled_system.Quantum_evolution(state_energy, time_list, survival_prob_list, electronic_survival_prob_list); // creat initial state (or read from file). Then complete simulation.

            state_energy_list.push_back(state_energy);
            time_list_all_states.push_back(time_list);
            survival_prob_list_all_states.push_back(survival_prob_list);
            electronic_survival_prob_list_all_states.push_back(electronic_survival_prob_list);
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

