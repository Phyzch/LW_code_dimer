//
// Created by phyzch on 6/23/20.
// Wrap up function to read data from input.txt and do output.
//
#include"../util.h"
#include"../system.h"
double nonadiabatic_coupling;

void full_system:: read_input_with_MPI(){
    if(my_id==0) {
        input.open(path + "input.txt"); // information recorded in input.txt
        if (!input.is_open()) {
            cout << "THE INFILE FAILS TO OPEN!" << endl;
            log<< "THE INFILE FAILS TO OPEN!" << endl;
            MPI_Abort(MPI_COMM_WORLD, -2);  // Abort the process and return error code -2. (input file can't open)
        }
        input >>  energy_window_size >> Rmax >> d.V_intra >> d.a_intra >> d.detector_energy_window_size >>  nonadiabatic_coupling;
        // read time used for simulation.  delt: time step.  tmax: maximum time for simulation.   tprint: time step to print result.
        input >> delt >> tmax >> tprint;
        // check if input is valid
        log.open(path + "log.txt");  // log to record the error information

        if (! log.is_open()){
            cout << "THE LOG FILE FAILS TO OPEN" << endl;
            MPI_Abort(MPI_COMM_WORLD, -1);      // Abort the process and return error code -1.  (log file can't open)
        }
        log <<"Number of Process running:  "<< num_proc <<endl;
        output.open(path+"output.txt");

        if( ! output.is_open()){
            cout<<"OUTPUT FILE FAILS TO OPEN."<<endl;
            log<< "OUTPUT FILE FAILS TO OPEN."<<endl;
            MPI_Abort(MPI_COMM_WORLD,-3); // Abort the process and return error code -3 (output file can't open).
        }

        if (delt <= 0 ) {
            cout << "Wrong time variable." << endl;
            log << "Wrong time variable." << "The simulation is cancelled." << endl;
            log.close();
            input.close();
            output.close();
            MPI_Abort(MPI_COMM_WORLD,-4);  // Abort with error code -4: Wrong time variable.
        }
        output << delt <<  " " << tmax << " " << tprint << endl;
    }
    // Broadcast hyper parameter to all process.
    // function:  int MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm)
    MPI_Bcast(&energy_window_size,1, MPI_DOUBLE, 0 , MPI_COMM_WORLD);
    MPI_Bcast(&Rmax,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&d.V_intra,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&d.detector_energy_window_size,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&nonadiabatic_coupling, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Bcast delt tmax tprint to other process.
    MPI_Bcast(&delt, 1, MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&tmax,1 ,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&tprint,1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // used to rescale the matrix element amplitude.
    cf = 0.0299792458*delt * pi2;
}


