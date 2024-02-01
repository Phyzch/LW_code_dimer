//
// Created by phyzch on 7/22/20.
//
#include"../util.h"
#include"../system.h"

// read parameters for system matrix set up
void system::read_MPI(ifstream &input, ofstream &output, ofstream &log) {
    // read energy level exciton_state_energy , initial wavefunction x_exciton, y_exciton
    int my_id;
    MPI_Comm_rank(MPI_COMM_WORLD,&my_id);
    if(my_id==0) {
        // only process 0 will do I/O here.
        // read number of detectors. typically 1 or 2.
        input >> exciton_state_num;
        if (exciton_state_num != 1 && exciton_state_num != 2) {
            log << "TLM NUMBER NOT SUPPORTED" << endl;
            input.close();
            log.close();
            output.close();
            exit(-5 ); // exciton_state_num is not right.
        }
        output << "exciton_state_num " << exciton_state_num << " ";
    }

    MPI_Bcast(&exciton_state_num, 1, MPI_INT, 0, MPI_COMM_WORLD);

    tlmatsize = exciton_state_num;  // system wave function array size.

    allocate_space();

    initialize_energy_level(input,output);

    initialize_wavefunction(input,output);

    initialize_state_energy();

};

void system::allocate_space(){
    // x_exciton: real part of system wave function
    // y_exciton: image part of system wave function
    // exciton_state_energy: energy level of system's eigen state, also diagonal part of our Hamiltonian matrix
    x_exciton = new double[exciton_state_num];
    y_exciton = new double[exciton_state_num];
    exciton_state_energy = new double[exciton_state_num];
}

void system::initialize_energy_level(ifstream & input, ofstream & output){
    // initialize energy of state
    int i;
    int my_id;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    if(my_id==0) {
        for (i = 0; i < exciton_state_num; i++) {
            input >> exciton_state_energy [i];
            output << exciton_state_energy [i] << " ";
        }
        output << endl;
    }
    // broadcast exciton_state_energy to all other process. (Variable in class system is all very small, do not have to data decomposition).
    MPI_Bcast(exciton_state_energy, exciton_state_num, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void system::initialize_wavefunction(ifstream & input, ofstream & output){
    //initialize wavefunction of state & normalize it
    int i;
    double norm = 0;

    if(my_id==0) {
        for (i = 0; i < tlmatsize; i++) {
            // x_exciton: real part of wave function.
            // y_exciton: imag part of wave function
            input >> x_exciton[i] >> y_exciton [i];
            output << x_exciton[i] << " " << y_exciton[i] << endl;
            norm = norm + pow(x_exciton[i], 2) + pow(y_exciton[i], 2);
        }
        norm = 1 / sqrt(norm);
        for (i = 0; i < tlmatsize; i++) {
            x_exciton[i] = x_exciton[i] * norm;
            y_exciton[i] = y_exciton[i] * norm;
        }
    }
    MPI_Bcast(x_exciton, tlmatsize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(y_exciton, tlmatsize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void system::initialize_state_energy(){
    // initialize energy of system state.
    int i,j;


}