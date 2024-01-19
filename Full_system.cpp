#include"system.h"
#include"util.h"
using namespace std;
// advise: When you get lost in all these functions. Keep eye on fullsystem::fullsystem() and full_system::Quantum_dynamics_evolution. Because
//they are functions do most of works and call other functions here.

// initialization of parameters and do some pre-coupling set up
full_system::full_system(string path1 , vector<vector<int>> & initial_state_quantum_number) {

	path = path1;
    d.path = path;
    // read hyper-parameter and time step from input.txt
    read_input_with_MPI();
    // here we study dimer model, therefore, there are two monomers.
    monomer_number = 2;
    d.monomer_number = monomer_number;

	s.read_MPI(input, output, log);
	d.read_MPI(input, output, s.exciton_state_num, path);
    d.construct_initial_state_MPI( initial_state_quantum_number);

    compute_monomer_vib_state_basis_set_size_MPI();

    d.construct_monomer_Hamiltonian_MPI(input, output, log, monomer1_vib_state_energy_all_pc, monomer2_vib_state_energy_all_pc,
                                        monomer_qn_list0,
                                        monomer_qn_list1);

    d.monomer1_vib_state_energy_all_pc = monomer1_vib_state_energy_all_pc;
    d.monomer2_vib_state_energy_all_pc = monomer2_vib_state_energy_all_pc;

    construct_dimer_Hamiltonian_matrix_with_energy_window_MPI();

	if(my_id ==0){
        cout<<"Finish constructing Matrix"<<endl;
        output_calculation_size_info(); // check if all matrix's dimension is right.
	}

}

// Doing Quantum Simulation with SUR algorithm, parallelized version.
void full_system::Quantum_dynamics_evolution(double & state_energy_for_record, vector<double> & time_list, vector<double> & survival_probability_list, vector<double> & electronic_state_survival_probability_list ,
                                             vector<vector<double>> & monomer_vib_energy) {
    // state_mode_list : record vibrational qn for initial state
    // time_list: record time.  survival probability list: record survival probability.  electronic_state_survival_probability_list : record electronic survival probability

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

	// Now we construct our wavefunction /phi for our monomer and full_system.
    d.initialize_monomer_state_MPI(log); // initialize monomer lower bright state

    Initial_state_MPI(); // construct initial state of whole system according to monomer state and system state.

    // prepare variable for evolution.
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
    double electronic_state_survival_prob;
    double electronic_survival_prob_sum;

    vector<double> electronic_state_label_array;
    generate_label_for_electronic_survival_prob_calculation(electronic_state_label_array);


	for (k = initial_step; k <= steps; k++) {

		if (k % psteps == 0) {
            // output result.
		    // Normalize wave function.
		    Normalize_wave_function();

            // update wave function received from other process.
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

            // --------------- code for computing survival probability in one potential energy surface
            electronic_state_survival_prob = 0;
            for(i=0;i<matsize;i++){
                electronic_state_survival_prob = electronic_state_survival_prob + (pow(x[i], 2) + pow(y[i], 2)) * electronic_state_label_array[i];
            }
            MPI_Allreduce(&electronic_state_survival_prob, &electronic_survival_prob_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            electronic_state_survival_probability_list.push_back(electronic_survival_prob_sum);
            // -----------------  end for code computing electronic survival probability.

            // code for computing vibrational energy in each monomer
            monomer1_vib_energy_each_pc = 0;
            monomer2_vib_energy_each_pc = 0;
            monomer1_vib_energy = 0;
            monomer2_vib_energy = 0;
            for(i = 0;i < matsize; i++){
                monomer1_vib_energy_each_pc = monomer1_vib_energy_each_pc + ( pow(x[i],2) + pow(y[i], 2)) * monomer1_vib_state_energy_all_pc[ vibrational_state_index_list[0][i] ];
                monomer2_vib_energy_each_pc = monomer2_vib_energy_each_pc + ( pow(x[i],2) + pow(y[i], 2)) * monomer2_vib_state_energy_all_pc[ vibrational_state_index_list[1][i] ];
            }
            MPI_Allreduce(&monomer1_vib_energy_each_pc, &monomer1_vib_energy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&monomer2_vib_energy_each_pc, &monomer2_vib_energy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            monomer1_vib_energy_list.push_back(monomer1_vib_energy);
            monomer2_vib_energy_list.push_back(monomer2_vib_energy);
            // end code for computing vib energy in each monomer

        }

        // evolve wave function one time step
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
	input.close();
	log.close();
	output.close();
	resource_output.close();

}
