#include"system.h"
#include"util.h"
using namespace std;
// advise: When you get lost in all these functions. Keep eye on fullsystem::fullsystem() and full_system::Quantum_evolution. Because
//they are functions do most of works and call other functions here.

// initialization of parameters and do some pre-coupling set up
full_system::full_system(string path1 , vector<vector<int>> & initial_state_quantum_number) {

	path = path1;
    d.path = path;
    // read hyper parameter and time step from input.txt
    read_input_with_MPI();

	s.read_MPI(input, output, log);
	d.read_MPI(input, output, s.electronic_state_num, path);
    d.construct_initial_state_MPI( initial_state_quantum_number);

    compute_detector_matrix_size_MPI();

    d.construct_dmatrix_MPI(input, output, log, dmat_diagonal_global0, dmat_diagonal_global1, vmode0, vmode1);

    d.dmat_diagonal_global0 = dmat_diagonal_global0;
    d.dmat_diagonal_global1 = dmat_diagonal_global1;

    construct_fullmatrix_with_energy_window_MPI();

	if(my_id ==0){
        cout<<"Finish constructing Matrix"<<endl;
        dimension_check(); // check if all matrix's dimension is right.
	}
}

// Doing Quantum Simulation with SUR algorithm, parallelized version.
void full_system::Quantum_evolution( double & state_energy_for_record, vector<double> & time_list, vector<double> & survival_probability_list, vector<double> & electronic_survival_probability_list ,
                                     vector<vector<double>> & monomer_vib_energy,
                                     vector<vector<vector<double>>> & EVD_electronic0_list,
                                     vector<vector<vector<double>>> & EVD_electronic1_list
                                     ) {
    // state_mode_list : record vibrational qn for initial state
    // time_list: record time.  survival probability list: record survival probability.  electronic_survival_probability_list : record electronic survival probability

    // ---------------- prepare variable for prepare evolution and call prepare_evolution ------------------------------
    int i, j, k,m;
    int steps, psteps;
    int irow_index, icol_index;

    // compute initial state energy.
    double initial_state_energy;
    if(my_id == initial_dimer_state_pc_id){
        initial_state_energy = mat[initial_dimer_state_index];
    }
    MPI_Bcast(&initial_state_energy, 1, MPI_DOUBLE, initial_dimer_state_pc_id, MPI_COMM_WORLD);

    // record the state energy
    state_energy_for_record = initial_state_energy;

    // record vib energy for each monomer
    vector<double> monomer1_vib_energy_list;
    vector<double> monomer2_vib_energy_list;
    double monomer1_vib_energy_each_pc, monomer2_vib_energy_each_pc;
    double monomer1_vib_energy, monomer2_vib_energy;

    // -----------------------------------------------------------------------------------
    // convert all matrix elements for ps-wavenumber units
    for (i = 0; i < matnum; i++) {
        mat[i] = cf * mat[i];
    }

	// Now we construct our wavefunction /phi for our detector and full_system.
    d.initialize_detector_state_MPI(log); // initialize detector lower bright state

    Initial_state_MPI(); // construct initial state of whole system according to detector state and system state.

    // prepare varibale for evolution.
    prepare_evolution();

    shift_mat();

    t = 0;

	steps = tmax / delt + 1; // Total number of steps for simulation.0
	psteps = tprint / delt;  // number of steps for printing out resuilt


	clock_t start_time, end_time, duration;
	start_time = clock();
	int initial_step = t / delt;

    // for computing survival probability
    double survival_prob = 0;
    double electronic_survival_prob;
    double electronic_survival_prob_sum;

    vector<double> electronic_state_label_array;
    generate_label_for_electronic_survival_prob_calculation(electronic_state_label_array);

    // quantum probability for coupled state
    vector<int> dimer_coupling_state_index_list;
    vector<int> dimer_coupling_state_pc_id_list;
    vector<vector<vector<int>>> dimer_coupled_state_quantum_number_list;
    vector<double>  dimer_coupling_state_energy_list;
    int coupled_dimer_state_num;
    double coupled_dimer_state_probability;

    construct_locally_coupled_states_for_monitor_Pt(dimer_coupling_state_index_list, dimer_coupling_state_pc_id_list, dimer_coupled_state_quantum_number_list, dimer_coupling_state_energy_list);
    coupled_dimer_state_num = dimer_coupled_state_quantum_number_list.size();

    // output coupled states.
    ofstream coupled_state_quantum_prob_output;
    if(my_id == 0){
        coupled_state_quantum_prob_output.open( path + "coupled_state_survival_prob.txt" );
        coupled_state_quantum_prob_output << coupled_dimer_state_num << endl;
        // output coupled state energy
        for(i=0;i<coupled_dimer_state_num;i++){
            coupled_state_quantum_prob_output << dimer_coupling_state_energy_list[i] << " ";
        }
        coupled_state_quantum_prob_output << endl;
        for(i=0;i<coupled_dimer_state_num;i++){
            for(m=0;m<d.electronic_state_num;m++){
                for(k=0;k<d.nmodes[m];k++){
                    coupled_state_quantum_prob_output << dimer_coupled_state_quantum_number_list[i][m][k] << " ";
                }
                coupled_state_quantum_prob_output << "       ";
            }
            coupled_state_quantum_prob_output << endl;
        }
    }

    // for exciton vibrational density
    mapping_mode0 = 0;
    mapping_mode1 = 1;
    // size: [1 + nmax[mapping_mode0], 1 + nmax[mapping_mode1]]
    vector<vector<double>> EVD_electronic0; // exciton vibrational density on electronic surface 0
    vector<vector<double>> EVD_electronic1; // exciton vibrational density on electronic surface 1

	for (k = initial_step; k <= steps; k++) {

		if (k % psteps == 0) {
            // output result.
		    // Normalize wave function.
		    Normalize_wave_function();
		    update_x_y();

            // record time
            time_list.push_back(t);


            // ---------- code for computing survival probability --------
            if(my_id == initial_dimer_state_pc_id){
                survival_prob = pow(x[initial_dimer_state_index],2) + pow(y[initial_dimer_state_index] , 2) ;
            }
            MPI_Bcast(&survival_prob,1, MPI_DOUBLE, initial_dimer_state_pc_id, MPI_COMM_WORLD);
            // record survival probability
            survival_probability_list.push_back(survival_prob);
            // --------- end for code computing survival probability ------

            // code for computing electronic survival probability
            electronic_survival_prob = 0;
            for(i=0;i<matsize;i++){
                electronic_survival_prob = electronic_survival_prob + ( pow(x[i],2) + pow(y[i], 2)) * electronic_state_label_array[i];
            }
            MPI_Allreduce(&electronic_survival_prob, &electronic_survival_prob_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            electronic_survival_probability_list.push_back(electronic_survival_prob_sum);

            //  end for code computing electronic survival probability.

            // code for computing vibrational energy in each monomer
            monomer1_vib_energy_each_pc = 0;
            monomer2_vib_energy_each_pc = 0;
            monomer1_vib_energy = 0;
            monomer2_vib_energy = 0;
            for(i=0;i<matsize;i++){
                monomer1_vib_energy_each_pc = monomer1_vib_energy_each_pc +  ( pow(x[i],2) + pow(y[i], 2)) *  dmat_diagonal_global0[ dstate[0][i] ];
                monomer2_vib_energy_each_pc = monomer2_vib_energy_each_pc +  ( pow(x[i],2) + pow(y[i], 2)) * dmat_diagonal_global1[ dstate[1][i] ];
            }
            MPI_Allreduce(&monomer1_vib_energy_each_pc, &monomer1_vib_energy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&monomer2_vib_energy_each_pc, &monomer2_vib_energy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            monomer1_vib_energy_list.push_back(monomer1_vib_energy);
            monomer2_vib_energy_list.push_back(monomer2_vib_energy);
            // end code for computing vib energy in each monomer

            // code for computing energy of coupled vibrational states
            if(my_id == 0){
                coupled_state_quantum_prob_output << t << endl;
            }
            for(i=0; i < coupled_dimer_state_num; i++){
                if(my_id == dimer_coupling_state_pc_id_list[i]){
                    coupled_dimer_state_probability = pow( x[dimer_coupling_state_index_list[i]] , 2 ) + pow( y[dimer_coupling_state_index_list[i]] ,2) ;
                }
                MPI_Bcast( & coupled_dimer_state_probability, 1, MPI_DOUBLE, dimer_coupling_state_pc_id_list[i], MPI_COMM_WORLD );
                if(my_id == 0){
                    coupled_state_quantum_prob_output << coupled_dimer_state_probability << " ";
                }
            }
            if(my_id == 0){
                coupled_state_quantum_prob_output << endl;
            }

            // code for computing probability of coupled vibrational states (reduced density matrix)
            compute_exciton_vibrational_density(EVD_electronic0, 0);
            compute_exciton_vibrational_density(EVD_electronic1, 1);
            EVD_electronic0_list.push_back(EVD_electronic0);
            EVD_electronic1_list.push_back(EVD_electronic1);

        }

        evolve_wave_func_one_step();


		t = t + delt;
	}

	end_time = clock();
	duration = end_time - start_time;
	if(my_id == 0) {
        log << "The total run time for parallel computing is " << (double(duration) /CLOCKS_PER_SEC)/60 << " minutes  for simulation time  " << tmax << endl;
        cout << "The total run time for parallel computing is " << (double(duration)/CLOCKS_PER_SEC)/60 << " minutes  for simulation time  " << tmax << endl;
    }

    monomer_vib_energy.push_back(monomer1_vib_energy_list);
    monomer_vib_energy.push_back(monomer2_vib_energy_list);

    // -------------------------- free space

    // ------------------------------------------------------------------------------------
    if (my_id == 0){
        input.close();
        log.close();
        output.close();
        resource_output.close();
        coupled_state_quantum_prob_output.close();
    }



}
