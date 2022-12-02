#include"system.h"
#include"util.h"
using namespace std;
// advise: When you get lost in all these functions. Keep eye on fullsystem::fullsystem() and full_system::Quantum_evolution. Because
//they are functions do most of works and call other functions here.

// initialization of parameters and do some pre-coupling set up
full_system::full_system(string path1) {

	path = path1;
    d.path = path;
    // read hyper parameter and time step from input.txt
    read_input_with_MPI();

	s.read_MPI(input, output, log);
	d.read_MPI(input, output, s.electronic_state_num, path);
    d.construct_initial_state_MPI(input, output);

    compute_detector_matrix_size_MPI();

    d.construct_dmatrix_MPI(input, output, log,dmat0,dmat1,vmode0,vmode1);

    construct_fullmatrix_with_energy_window_MPI();

	if(my_id ==0){
        cout<<"Finish constructing Matrix"<<endl;
        dimension_check(); // check if all matrix's dimension is right.
	}
}

// Doing Quantum Simulation with SUR algorithm, parallelized version.
void full_system::Quantum_evolution() {
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

    // vibrational survival prob
    ofstream survival_prob_out;
    if(my_id == 0){
        survival_prob_out.open(path + "survival_prob.txt");
        // output energy

        survival_prob_out << initial_state_energy << " "; // energy for state
        survival_prob_out << endl;

        // output mode quanta
        for(m=0;m<d.electronic_state_num;m++){
            for(j=0;j<d.nmodes[m];j++){
                survival_prob_out << d.initial_detector_state[m][j]<<"  ";
            }
            survival_prob_out << endl;
        }
    }

    // electronic survival probability
    ofstream electronic_survival_prob_out;
    if(my_id == 0){
        electronic_survival_prob_out.open(path + "electronic_survival_prob.txt");
        // output energy

        electronic_survival_prob_out << initial_state_energy << endl;

        // output mode quanta
        for(m=0;m<d.electronic_state_num;m++){
            for(j=0;j<d.nmodes[m];j++){
                electronic_survival_prob_out << d.initial_detector_state[m][j]<<"  ";
            }
            electronic_survival_prob_out << endl;
        }
    }


	for (k = initial_step; k <= steps; k++) {

		if (k % psteps == 0) {
            // output result.
		    // Normalize wave function.
		    Normalize_wave_function();
		    update_x_y();

            // ---------- code for computing survival probability --------
            if(my_id == initial_dimer_state_pc_id){
                survival_prob = pow(x[initial_dimer_state_index],2) + pow(y[initial_dimer_state_index] , 2) ;
            }
            MPI_Bcast(&survival_prob,1, MPI_DOUBLE, initial_dimer_state_pc_id, MPI_COMM_WORLD);
            if(my_id == 0){
                // output survival prob
                survival_prob_out << t << endl;
                survival_prob_out << survival_prob << endl;
            }
            // --------- end for code computing survival probability ------

            // code for computing electronic survival probability
            electronic_survival_prob = 0;
            for(i=0;i<matsize;i++){
                electronic_survival_prob = electronic_survival_prob + ( pow(x[i],2) + pow(y[i], 2)) * electronic_state_label_array[i];
            }
            MPI_Allreduce(&electronic_survival_prob, &electronic_survival_prob_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            if(my_id == 0){
                // output electronic survival probability
                electronic_survival_prob_out << t << endl;
                electronic_survival_prob_out << electronic_survival_prob_sum << endl;
            }

            //  end for code computing electronic survival probability.
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


    // -------------------------- free space

    // ------------------------------------------------------------------------------------
	input.close();
	log.close();
	output.close();
	resource_output.close();
    survival_prob_out.close();
    electronic_survival_prob_out.close();
}
