//
// Created by phyzch on 7/22/20.
//
#include"../util.h"
#include"../system.h"

// read parameters for system matrix set up
void system::read_MPI(ifstream &input, ofstream &output, ofstream &log) {
    // read energy level electronic_state_energy , initial wavefunction x_electronic, y_electronic
    int my_id;
    MPI_Comm_rank(MPI_COMM_WORLD,&my_id);
    if(my_id==0) {
        // only process 0 will do I/O here.
        // read number of detectors. typically 1 or 2.
        input >> electronic_state_num;
        if (electronic_state_num != 1 && electronic_state_num != 2) {
            log << "TLM NUMBER NOT SUPPORTED" << endl;
            input.close();
            log.close();
            output.close();
            exit(-5 ); // electronic_state_num is not right.
        }
        output << "electronic_state_num " << electronic_state_num << " ";
    }
    MPI_Bcast(&electronic_state_num, 1, MPI_INT, 0, MPI_COMM_WORLD);
    tlmatsize = electronic_state_num;  // system wave function array size.

    initialize_energy_level(input,output);

    initialize_wavefunction(input,output);

    initialize_state_energy();

};


void system::initialize_energy_level(ifstream & input, ofstream & output){
    // initialize energy of state
    int i;
    int my_id;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    if(my_id==0) {
        for (i = 0; i < electronic_state_num; i++) {
            input >> electronic_state_energy [i];
            output << electronic_state_energy [i] << " ";
        }
        output << endl;
    }
    // broadcast electronic_state_energy to all other process. (Variable in class system is all very small, do not have to data decomposition).
    MPI_Bcast(electronic_state_energy, electronic_state_num, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void system::initialize_wavefunction(ifstream & input, ofstream & output){
    //initialize wavefunction of state & normalize it
    int i;
    double norm = 0;

    if(my_id==0) {
        for (i = 0; i < tlmatsize; i++) {
            input >> x_electronic[i] >> y_electronic[i];
            output << x_electronic[i] << " " << y_electronic[i] << endl;
            norm = norm + pow(x_electronic[i], 2) + pow(y_electronic[i], 2);
        }
        norm = 1 / sqrt(norm);
        for (i = 0; i < tlmatsize; i++) {
            x_electronic[i] = x_electronic[i] * norm;
            y_electronic[i] = y_electronic[i] * norm;
        }
    }
    MPI_Bcast(x_electronic, tlmatsize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(y_electronic, tlmatsize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void system::initialize_state_energy(){
    // initialize energy of system state.
    int i,j;

    for (i = 0; i < tlmatsize; i++) {
        tlmat[i] = electronic_state_energy[i];
    }

}