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
	friend class detector;
	double *x_electronic, *y_electronic, *electronic_state_energy, *tlmat; //
	int electronic_state_num, tlmatsize;
	void read_MPI(ifstream &input, ofstream &output, ofstream &log);
    void initialize_energy_level(ifstream & input, ofstream & output);
    void initialize_wavefunction(ifstream & input, ofstream & output);
    void initialize_state_energy();
    void allocate_space();
	system();
	~system();
};

class detector {
private:
    const int fillfrac = 10;
    const static int detdim =80;  // maximum n*n dimension of our detector matrix
    int dmatdim; // detector matrix size
    int electronic_state_num;
    string path;

    int * initial_state_index;
    int * initial_state_pc_id;
public:
	// for mode:  modtype: =0 dark or =1 bright state.  nmax: maximum state number for each mode. nmodes: total mode numbers.
	friend class full_system;
	friend class system;
	int *nmodes, **nmax;
	int *dmatsize;
    int *dmatnum , *doffnum;  // detector matrix elemetn array
	vector<int> total_dmat_size; // size of whole matrix across various process.
	vector <int> total_dmat_num, total_dmat_off_num; // total matrix number and off-diagonal number across various process.
    vector<vector <int>> * dv_all;
	int ** dmatsize_each_process;
	int ** doffnum_each_process;
	int ** dmatnum_each_process;  // record detector matrix element number in each process.
	int ** dmat_offset_each_process; // record local first detector matrix's index in global matrix.

	double ** total_dmat;
	int ** total_dirow, ** total_dicol; // dirow, dicol, dmat in all process.

	vector<vector<int>> *dv;  //dv: the q.n. for states (m,i) at coordinate j.  [m][i][j]: [detector][state][mode]
	vector<int> *dirow;
	vector<int> *dicol;
	int *deln;  // deln= |n_{i} - n_{j}| at coordinate k
	double *nbar;
	double **mfreq; // frequency of mode
	double **aij;
	vector <double> * xd, * yd; // wavefunction of detector state
    double ** xd_all, ** yd_all;
	vector<double> * dmat; // matrix
	double *proptime; // we set proptime for two different mode to be the same

    // for EV coupling
    double ** electron_phonon_coupling;
    double **** franck_condon_factor_table;

    int ** remoteVecCount,  ** remoteVecPtr,  **  remoteVecIndex;
    int ** tosendVecCount,  **tosendVecPtr,  ** tosendVecIndex;
    int * to_send_buffer_len, * to_recv_buffer_len;
    double ** recv_xd,  ** recv_yd,  ** send_xd,  ** send_yd;
    vector <int > * local_dirow;
    vector <int > * local_dicol;


    //matflag	: = 1 compute detector matrixusing Madsen scaling Hamiltonian
    //maxdis : largest distance in q.n.space for which matrix element is calculated
    //cutoff : perturbation cutoff criterion V / delta - E(typ. 0.05) for states within a single detector
    //cutoff2 : same for states between different detectors : a(i)a(j) / delta - E must be greater than cutoff2
    //kelvin : detector temperature in kelvin
	int   maxdis;
	double cutoff;
    double Franck_condon_factor_cutoff;
    vector<vector<vector<int>>> nonadiabatic_coupled_d_state; // state nonadiabatically coupled to given d state.
    vector<vector<vector<double>>> nonadiabatic_coupled_d_state_franck_condon; // franck_condon factor for nonadiabtatically coupled state

	double V_intra, a_intra; // intra detector coupling strength.  a_intra = coupling strength for mfreq = 50.
    double detector_energy_window_size;
	int  ** initial_detector_state; // record bright mode for two detectors when we try to see decoherence in our model.
	double * initial_state_energy;
	detector();
	~detector();
	void allocate_space();
    void allocate_space_single_detector(int detector_index);
	void read_MPI(ifstream & input, ofstream & output, int electronic_state_num1, string path);


    // MPI version of function.
    void construct_dv_dirow_dicol_dmatrix_MPI(ofstream & log, vector<double> & dmat0,  vector<double> & dmat1,  vector<vector<int>> & vmode0, vector<vector<int>> & vmode1);
    void construct_dmatrix_MPI(ifstream & input, ofstream & output, ofstream & log, vector<double> & dmat0,  vector<double> & dmat1,  vector<vector<int>> & vmode0, vector<vector<int>> & vmode1);
    void compute_detector_offdiag_part_MPI(ofstream & log,vector<double> & dmat0,  vector<double> & dmat1,  vector<vector<int>> & vmode0, vector<vector<int>> & vmode1);
    void broadcast_total_dmat();  // broadcast dmat, dirow , dicol to form total_dirow, total_dicol, total_dmat
    void broadcast_dmatnum_doffnum(); // broadcast dmatnum, doffnum, dmatnum_each_process, doffnum_each_process, total_dmatnum etc.
    void gather_xd_yd();
    void Scatter_xd_yd();
    void Scatter_dv(vector<int> & total_mat_num);
    void Scatter_dirow_dicol_dmatrix(vector <double> * dmat_all, vector<int> * dirow_data, vector<int> * dicol_data, int ** vector_size,int **vector_displs , ofstream & log);
    int construct_receive_buffer_index(int * remoteVecCount_element, int * remoteVecPtr_element, int * remoteVecIndex_element, int detector_index);
    void prepare_evolution();
    // MPI version of SUR for one detector for each timestep.
    void update_dx_dy(int detector_index);
    void SUR_onestep_MPI( int detector_index, double cf);
    void construct_initial_state_MPI(ifstream & input, ofstream & output);
    void initialize_detector_state_MPI(ofstream & log);

    // used to broadcast dv_all , vmode0, vmode1 , dmat0, dmat1
    void Broadcast_dv_all();

    void compute_important_state_index();

    // for nonadiabatic franck condon factor
    void compute_franck_condon_factor_table();
    void find_franck_condon_factor_for_monomer_states();

};

class full_system {
	// detector+ system
private:
	int matdim;  // matdim is maximum size of full matrix
	int matnum, offnum, matsize; // matnum is total matrix element number, it should be smaller than matdim.
	                                // in our program, usually we set matdim=matnum and use these variable interchangably.
								 //offnum: off diagonal matrix element number. matsize: matrix size for full matrix (number of diagonal element)
    int total_matnum, total_offnum, total_matsize;
    int * matsize_each_process, *mat_offnum_each_process, *matnum_each_process;
    int * matsize_offset_each_process, * matnum_offset_each_process;
    vector <double> x;
    vector <double> y; // x[matsize], y[matsize]
	vector <double> mat; // full system matrix, size: matdim
	vector <int> irow, icol; // row index and column index of matrices, size:matdim
	vector <int>sstate;
	vector<int> * dstate;  // sstate[matsize], dstate[matsize].  index of system state and detector for that matrix element. (This could be replaced by function to save space)

	double total_energy;
	double norm; // used to check normalization
    double total_norm;
				 // below are some variables for density matrix

    vector <quotient_state> d1list;  // state in quotient Hilbert space for detector 1
    vector <quotient_state> d2list;  // state in quotient Hilbert space for detector 2

    vector <sys_quotient_state> slist;  // state in quotient Hilbert space for system.

    // vmode,dmat for each detector.
    vector<vector<int>> vmode0;
    vector<vector<int>> vmode1;
    vector<double> dmat0;
    vector<double> dmat1;


    // used for receiving and sending vedtor x , y from/to other process
    int * remoteVecCount,  *remoteVecPtr,  * remoteVecIndex,
     * tosendVecCount, * tosendVecPtr,  * tosendVecIndex;
    double * recv_x,  * recv_y, * send_x, * send_y;
    int  to_send_buffer_len,  to_recv_buffer_len;
    vector<int>  local_irow;
    vector <int>  local_icol;

    // used for recving and sending vector x,y for computing system_energy

public:
	class system s;
	class detector d;

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

	full_system(string path1);
	~full_system();
	void dimension_check();
	void Quantum_evolution();;

    // MPI version of code:
    void read_input_with_MPI();
    void compute_detector_matrix_size_MPI();
    void pre_coupling_evolution_MPI(int initial_state_choice);
    void construct_fullmatrix_with_energy_window_MPI();
    void compute_sstate_dstate_diagpart_dirow_dicol_MPI( );
    void construct_quotient_state_all_MPI();
    vector<vector<quotient_state>>  Gather_quotient_state_vec();
    vector<quotient_state> sort_d1list(vector<vector<quotient_state>> & d1list_each_proces);
    void Scatter_sorted_d1list(vector<quotient_state> & sorted_d1list);
    void construct_q_index_MPI();
    void rearrange_d1list();
    void compute_offdiagonal_part_MPI();
    void compute_dmat_off_diagonal_matrix_in_full_matrix_one_dmat_MPI(int index,vector < double > & mat,
            vector<int> & irow, vector<int> & icol);
    void compute_dmat_off_diagonal_matrix_in_full_matrix_MPI(vector < double > & mat,vector  <int> & irow, vector<int> & icol);

    void rearrange_off_diagonal_term(vector < double > & mat,vector  <int> & irow, vector<int> & icol);
    void  combine_offdiagonal_term(vector<double> & d_off_mat, vector<int> & d_off_irow, vector<int> & d_off_icol);

    void Initial_state_MPI();

    // function to prepare and evolve photon+detector system:
    int construct_recvbuffer_index();
    void prepare_evolution();
    void  update_x_y();
    void update_x();
    void update_y();
    void evolve_wave_func_one_step();
    void full_system_SUR_one_step();
    // Output function MPI version

    void etot_MPI(double * hx, double * hy);

    void shift_mat();

    void gather_x_y(double * x_all, double * y_all); // gather x,y to save the state or load the state to file.
    void scatter_x_y(double * x_all, double * y_all); // scatter x_all, y_all to x,y.

    void Normalize_wave_function();

    // search x_index using all index in electronic dof, index in monomer1 and index in monomer2
    int search_full_sys_matrix_given_sd_matrix(int s_state, int d1_state, int d2_state); // find index in full matrix given index in electronic_state, in monomer1, and in monomer2.


};



