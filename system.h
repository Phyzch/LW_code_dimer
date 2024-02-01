#pragma once
#include<iostream>
#include<iomanip>
#include <fstream>
#include <complex>
#include <vector>
#include<ctime>
#include"quotient_state.h"
using namespace std;

// Information for MPI program
extern int my_id;
extern int num_proc;

class system {
public:
	friend class full_system;
	friend class monomer;
	double *x_exciton, *y_exciton, *exciton_state_energy; //
	int exciton_state_num, tlmatsize;
	void read_MPI(ifstream &input, ofstream &output, ofstream &log);
    void initialize_energy_level(ifstream & input, ofstream & output);
    void initialize_wavefunction(ifstream & input, ofstream & output);
    void initialize_state_energy();
    void allocate_space();
	system();
	~system();
};

class monomer {
private:
    int exciton_state_num;

    string path;

    int * initial_state_index;
    int * initial_state_pc_id;
public:
    int monomer_number;
	// for mode:  modtype: =0 dark or =1 bright state.  nmax: maximum state number for each mode. nmodes: total mode numbers.
	friend class full_system;
	friend class system;
	int *nmodes, **nmax;
	int *monomer_matsize;
    int *monomer_matnum , *monomer_offnum;  // monomer matrix elemetn array
	vector<int> total_monomer_mat_size; // size of whole matrix across various process.
	vector <int> total_monomer_mat_num, total_monomer_mat_off_num; // total matrix number and off-diagonal number across various process.
    vector<vector <int>> * monomer_vibrational_states_all;
	int ** monomer_matsize_each_process;
	int ** monomer_offnum_each_process;
	int ** monomer_matnum_each_process;  // record monomer matrix element number in each process.
    int ** monomer_matsize_offset_each_process;
	int ** monomer_mat_offset_each_process; // record local first monomer matrix's index in global matrix.

	double ** total_monomer_mat;
	int ** total_monomer_irow, ** total_monomer_icol; // monomer_irow, monomer_icol, monomer_mat in all process.
    vector<double> monomer1_vib_state_energy_all_pc;
    vector<double> monomer2_vib_state_energy_all_pc;

	vector<vector<int>> *monomer_vibrational_states_quantum_number_list;  //monomer_vibrational_states_quantum_number_list: the q.n. for states (m,i) at coordinate j.  [m][i][j]: [monomer][state][mode]

    vector<int> *monomer_irow;
	vector<int> *monomer_icol;
    vector<double> * monomer_mat; // matrix

	int *deln;  // deln= |n_{i} - n_{j}| at coordinate k
	double *nbar;
	double **mfreq; // frequency of mode
	double **aij;
	vector <double> * xd, * yd; // wavefunction of monomer state
    double ** xd_all, ** yd_all;

	double *proptime; // we set proptime for two different mode to be the same

    // for EV coupling
    double ** exciton_phonon_coupling;
    double **** franck_condon_factor_table;

    int ** remoteVecCount,  ** remoteVecPtr,  **  remoteVecIndex;
    int ** tosendVecCount,  **tosendVecPtr,  ** tosendVecIndex;
    int * to_send_buffer_len, * to_recv_buffer_len;
    double ** recv_xd,  ** recv_yd,  ** send_xd,  ** send_yd;
    vector <int > * local_dirow;
    vector <int > * local_dicol;


    //matflag	: = 1 compute monomer matrixusing Madsen scaling Hamiltonian
    //maxdis : largest distance in q.n.space for which matrix element is calculated
    //cutoff : perturbation cutoff criterion V / delta - E(typ. 0.05) for states within a single monomer
    //cutoff2 : same for states between different detectors : a(i)a(j) / delta - E must be greater than cutoff2
    //kelvin : monomer temperature in kelvin
	int   maxdis;
	double cutoff;
    double Franck_condon_factor_cutoff;
    vector<vector<vector<int>>> nonadiabatic_coupled_monomer_state; // state nonadiabatically coupled to given d state.
    vector<vector<vector<double>>> nonadiabatic_coupled_monomer_state_franck_condon; // franck_condon factor for nonadiabtatically coupled state

	double V_intra, a_intra; // intra monomer coupling strength.  a_intra = coupling strength for mfreq = 50.
    double vibrational_energy_window_size;
	int  ** initial_vibrational_state; // record bright mode for two detectors when we try to see decoherence in our model.
	double * initial_state_energy;
	monomer();
	~monomer();
	void allocate_space();
    void allocate_space_single_detector(int detector_index);
	void read_MPI(ifstream & input, ofstream & output, int electronic_state_num1, string path);


    // MPI version of function.
    void construct_dv_dirow_dicol_dmatrix_MPI(ofstream & log, vector<double> & dmat0,  vector<double> & dmat1,  vector<vector<int>> & vmode0, vector<vector<int>> & vmode1);
    void construct_monomer_Hamiltonian_MPI(ifstream & input, ofstream & output, ofstream & log, vector<double> & dmat_diagonal_global0, vector<double> & dmat_diagonal_global1, vector<vector<int>> & vmode0, vector<vector<int>> & vmode1);
    void compute_monomer_offdiag_part_MPI(ofstream & log, vector<double> & dmat0, vector<double> & dmat1, vector<vector<int>> & vmode0, vector<vector<int>> & vmode1);
    void broadcast_total_dmat();  // broadcast monomer_mat, monomer_irow , monomer_icol to form total_monomer_irow, total_monomer_icol, total_monomer_mat
    void broadcast_dmatnum_doffnum(); // broadcast monomer_matnum, monomer_offnum, monomer_matnum_each_process, monomer_offnum_each_process, total_dmatnum etc.
    void gather_xd_yd();
    void Scatter_xd_yd();
    void Scatter_dv(vector<int> & total_mat_num);
    void Scatter_dirow_dicol_dmatrix(vector <double> * dmat_all, vector<int> * dirow_data, vector<int> * dicol_data, int ** vector_size,int **vector_displs , ofstream & log);
    int construct_receive_buffer_index(int * remoteVecCount_element, int * remoteVecPtr_element, int * remoteVecIndex_element, int detector_index);
    void prepare_evolution();
    // MPI version of SUR for one monomer for each timestep.
    void update_dx_dy(int detector_index);
    void SUR_onestep_MPI( int detector_index, double cf);
    void construct_initial_state_MPI(vector<vector<int>> & initial_state_quantum_num);
    void initialize_monomer_state_MPI(ofstream & log);

    // used to broadcast monomer_vibrational_states_all , monomer_qn_list0, monomer_qn_list1 , monomer1_vib_state_energy_all_pc, monomer2_vib_state_energy_all_pc
    void Broadcast_dv_all();

    void compute_initial_vibrational_state_index();

    // for nonadiabatic franck condon factor
    void compute_franck_condon_factor_table();
    void find_franck_condon_factor_for_monomer_states();

    // compute N loc for each monomer
    void compute_local_density_of_state_subroutine(int monomer_index, vector<double> & state_energy_global_matrix,
                                                   vector<int> & coupling_state_index_list, vector<vector<int>> & coupling_state_qn_list,
                                                   vector<double> & coupling_state_strength, vector<double> & coupling_state_energy_diff,
                                                   double & effective_coupling_number);

    void compute_local_density_of_state(vector<vector<vector<int>>> & coupling_state_index_list,
                                        vector<vector<vector<vector<int>>>> & coupling_state_qn_list,
                                        vector<vector<vector<double>>> & coupling_state_strength_list,
                                        vector<vector<vector<double>>> & coupling_state_energy_diff_list,
                                        vector<vector<double>> & effective_coupling_number_list);

};

class full_system {
	// monomer+ system
private:
    int monomer_number;
	int matnum, offnum, matsize; // matnum is total matrix element number, it should be smaller than matdim.
	                                // in our program, usually we set matdim=matnum and use these variable interchangably.
								 //offnum: off diagonal matrix element number. matsize: matrix size for full matrix (number of diagonal element)
    int total_matnum, total_offnum, total_matsize;
    int * matsize_each_process, *mat_offnum_each_process, *matnum_each_process;
    int * matsize_offset_each_process, * matnum_offset_each_process;
    vector <double> real_part_wave_func;
    vector <double> imag_part_wave_func; // real_part_wave_func[matsize], imag_part_wave_func[matsize]
	vector <double> mat; // full system matrix, size: matdim
	vector <int> irow, icol; // row index and column index of matrices, size:matdim
	vector <int> exciton_state_index_list;
	vector<int> * vibrational_state_index_list;  // exciton_state_index_list[matsize], vibrational_state_index_list[matsize].  index of system state and monomer for that matrix element. (This could be replaced by function to save space)

    vector<int> exciton_state_index_list_all;
    vector<int> * vibrational_state_index_list_all;

    int initial_dimer_state_index; // index in process (local index)
    int initial_dimer_state_pc_id; // process id for index.

	double total_energy;
	double norm; // used to check normalization
    double total_norm;
				 // below are some variables for density matrix

    vector <quotient_state> monomer1_quotient_state_list;  // state in quotient Hilbert space for monomer 1
    vector <quotient_state> monomer2_quotient_state_list;  // state in quotient Hilbert space for monomer 2


    // vmode,monomer_mat for each monomer.
    vector<vector<int>> monomer_qn_list0;
    vector<vector<int>> monomer_qn_list1;
    vector<double> monomer1_vib_state_energy_all_pc;
    vector<double> monomer2_vib_state_energy_all_pc;


    // used for receiving and sending vedtor real_part_wave_func , imag_part_wave_func from/to other process
    int * remoteVecCount,  *remoteVecPtr,  * remoteVecIndex,
     * tosendVecCount, * tosendVecPtr,  * tosendVecIndex;
    double * recv_real_wave_func,  * recv_imag_wave_func, * send_real_wave_func, * send_imag_wave_func;
    int  to_send_buffer_len,  to_recv_buffer_len;
    vector<int>  local_irow;
    vector <int>  local_icol;

    // used for recving and sending vector real_part_wave_func,imag_part_wave_func for computing system_energy

public:
	class system s;
	class monomer d;

	// output and input of file
	ofstream output; // output result we record.
	ofstream log;  // output status (problem with code)
	ofstream resource_output;
	ifstream input;
	ofstream save;
	ifstream load;
	// timestep variable
	double delt, tmax, tprint;
	double t; // check the time

	double cf; // energy scale

	string path;

	full_system(string path1 , vector<vector<int>> & initial_state_quantum_number);
	~full_system();
	void output_calculation_size_info();
	void Quantum_dynamics_evolution(double & state_energy_for_record, vector<double> & time_list, vector<double> & survival_probability_list, vector<double> & electronic_state_survival_probability_list ,
                                    vector<vector<double>> & monomer_vib_energy);;

    // MPI version of code:
    void read_input_with_MPI();
    void compute_monomer_vib_state_basis_set_size_MPI();
    void pre_coupling_evolution_MPI(int initial_state_choice);
    void construct_dimer_Hamiltonian_matrix_with_energy_window_MPI();
    void compute_sstate_dstate_diagpart_dirow_dicol_MPI( );
    void construct_quotient_state_all_MPI();
    vector<vector<quotient_state>>  Gather_quotient_state_vec();
    vector<quotient_state> sort_d1list(vector<vector<quotient_state>> & d1list_each_proces);
    void Scatter_sorted_d1list(vector<quotient_state> & sorted_d1list);
    void construct_anharmonic_coupling_info_index_list_MPI();
    void rearrange_monomer1list();
    void compute_full_Hamiltonian_offdiagonal_part_MPI();
    void compute_anharmonic_coupling_in_full_matrix_in_one_monomer_MPI(int monomer_index, vector < double > & anharmonic_coupling_mat,
                                                                       vector<int> & anharmonic_coupling_irow, vector<int> & anharmonic_coupling_icol);
    void compute_monomer_anharmonic_coupling_in_full_matrix_MPI(vector < double > & anharmonic_coupling_mat, vector  <int> & anharmonic_coupling_irow, vector<int> & anharmonic_coupling_icol);

    void compute_nonadiabatic_offdiagonal_matrix_full_system(vector < double > & nonadiabatic_off_mat,
                                                             vector  <int> & nonadiabatic_off_irow,
                                                             vector<int> & nonadiabatic_off_icol);

    void rearrange_matrix_element_in_different_pc(vector < double > & mat, vector  <int> & irow, vector<int> & icol);
    void  combine_offdiagonal_term(vector<double> & anharmonic_coupling_mat, vector<int> & anharmonic_coupling_irow, vector<int> & anharmonic_coupling_icol,
                                   vector<double> & nonadiabatic_off_mat, vector<int> & nonadiabatic_off_irow, vector<int> & nonadiabatic_off_icol);

    void Initial_state_MPI();

    // function to prepare and evolve photon+monomer system:
    int construct_recvbuffer_index();
    void prepare_evolution();
    void  update_x_y();
    void update_real_part();
    void update_imag_part();
    void evolve_wave_func_one_step();
    void full_system_SUR_one_step();
    // Output function MPI version

    void etot_MPI(double * hx, double * hy);

    void shift_mat();

    void gather_x_y(double * x_all, double * y_all); // gather real_part_wave_func,imag_part_wave_func to save the state or load the state to file.
    void scatter_x_y(double * x_all, double * y_all); // scatter x_all, y_all to real_part_wave_func,imag_part_wave_func.

    void Normalize_wave_function();

    // search x_index using all index in electronic dof, index in monomer1 and index in monomer2
    int search_full_hamiltonian_state_index(int exciton_state, int monomer1_state, int monomer2_state); // find index in full matrix given index in electronic_state, in monomer1, and in monomer2.

    // for computing electronic survival probability
    void generate_label_for_electronic_survival_prob_calculation(vector<double> & electronic_state_label_array);

    // compute Nloc for each monomer.
    void compute_local_density_of_state(vector<vector<vector<int>>> & coupling_state_index_list,
                                        vector<vector<vector<vector<int>>>> & coupling_state_qn_list,
                                        vector<vector<vector<double>>> & coupling_state_strength_list,
                                        vector<vector<vector<double>>> & coupling_state_energy_diff_list,
                                        vector<vector<double>> & effective_coupling_number_list);
};



