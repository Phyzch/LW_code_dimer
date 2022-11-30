#include"system.h"
#include"util.h"
using namespace std;
// advise: When you get lost in all these functions. Keep eye on fullsystem::fullsystem() and full_system::Quantum_evolution. Because
//they are functions do most of works and call other functions here.

double energy_window_size;  // size of energy window, we will only include whole system state whose energy difference with

double detector_coupling_time = 20 ; // time for detector couple to each other. Beyond that time, detector_coupling strength will be 0.

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


	// Here I'd like to create a new file to output detector reduced density matrix and average quanta in each mode.  You can comment this code if you don't want this one.
	if(my_id==0) {
        Detector_output.open(path + "Detector_output.txt"); // output the information for next simulation.
        Detector_mode_quanta.open(path + "Detector_mode_quanta.txt");

        log << "Start SUR Calculation" << endl;
        log << "Total Time steps:" << steps << endl;
        output<< "time    s1    s2      Trsr2    se      de[0]       de[1]"
                 "    e    norm   s0 "<< endl;
    }


	clock_t start_time, end_time, duration;
	start_time = clock();
	int initial_step = t / delt;


	for (k = initial_step; k <= steps; k++) {
		if (k % psteps == 0) {
		    Normalize_wave_function();
		    update_x_y();
			if(my_id==0) {
                output << "Steps: " << k << endl;
            }
        }

		t = t + delt;

		// update tosend_xd for sending  to other process
        update_y();
        // SUR algorithm
        for(i=0;i<matnum;i++){
            irow_index = local_irow[i];
            icol_index= local_icol[i];
            x[irow_index] = x[irow_index] + mat[i] * y[icol_index];
        }
        update_x();
        for(i=0;i<matnum;i++){
            irow_index=local_irow[i];
            icol_index= local_icol[i];
            y[irow_index] = y[irow_index] - mat[i] * x[icol_index];
        }
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
	Detector_output.close();
	resource_output.close();
}
