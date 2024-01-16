# Note about code

## How to configure and run the code:

1. Make sure folder name of project is the same as that appear in CMAKEList.txt.
   
   For example if in CMAKEList.txt you write:
   
   (1) project(  my_project_LW_model )

        Then the folder name must be : **my_project_LW_model**        

            



2. in util.h

in my machine, I have #include<mpi/mpi.h>. This may be problematic on cluster or your computer.  May change it to #include <mpi.h>

3. ```bash
   cd SOURCE_DIRECTORY/Release 
   cmake ..
   make
   mpirun -np 10 ./executable (10 is number of process you want to use
   change it to different number if you want to use more or less)
   ```
   
   4. How to compile the project
      
      check : [c++ - Debug vs Release in CMake - Stack Overflow](https://stackoverflow.com/questions/7724569/debug-vs-release-in-cmake)
      
      ```bash
      # create file for debugging
      mkdir Debug
      cd ./Debug
      cmake -DCMAKE_BUILD_TYPE=Debug ..
      make 
      
      # create file for computing (production usage)
      mkdir Release
      cd Release
      cmake -DCMAKE_BUILD_TYPE=Release ..
      make
      ```
   
   ## Structure of Code
   
   We follow the order that program execute the functions :
   
   ### Full_system(string path1)

```cpp
full_system::full_system(string path1) {

 path = path1;
 d.path = path;
 d.cvpt_path = cvpt_path1;

 // read parameter and time step from input.txt
 read_input_with_MPI();

 s.read_MPI(input, output, log);
 d.read_MPI(input, output, log, s.tlnum, s.tldim,path);
 d.construct_initial_state_MPI(input,output);

 // We construct vibrational state for each electronic state . vibrational state is represented by their mode quanta , for example (1,0,0,0,2,0,1) etc.
 compute_detector_matrix_size_MPI_sphere();

 // we construct anharmonic coupling between states, thus finish constructing Hamiltonian.
 d.construct_monomer_Hamiltonian_MPI(input,output,log,dmat0,dmat1,monomer_qn_list0,monomer_qn_list1);
 cout<<"Finish constructing Matrix"<<endl;

}
```

This function construct Hamiltonian  $H$ for our systems.

#### 1 . s.read_MPI()  , d.read_MPI()

s.read_MPI()  , d.read_MPI()  : read parameter from input file 

(only process with rank 0 read the parameters and then broadcast parameter to all other process.  

- **You can't let multiple process access one file at same time, this may cause error.**

- **Make sure parameter read from input.txt is broadcasted to all other process, otherwise the simulation in other process will raise error because they don't know the value of the parameter !!!**

)

There are many paramters read but not used there, You can delete it as you wish.

#### 2. construct_initial_state_MPI()

construct initial state.

- First index [0] : stands for electronic state  (0 for electronic state 0, 
  
  1 for electronic state 1)

    |01111> : at electronci state 0, and vibrational state is |1111>

|10101> : at electronic state 1, vib state is |0101>

        We set $|n_{0} = 0 \rangle $ stands for ground electronic state, 

                   $|n_{0} = 1 \rangle $ stands for excited electronic state 

- Second index [1] : stands for first vibrational state. 

In our model, first vibrational state is shifted by $\lambda (b + b^{\dagger})$ , and $b \rightarrow b^{\dagger}$ as in Hamiltonian : 

$H^{\prime \prime}=\sum_{\sigma}|\sigma\rangle\left\langle\sigma\left|\left(\frac{\sigma}{2} \epsilon+\sigma \sum_{\alpha} \lambda_{\alpha \sigma}\left(b_{\alpha}^{\dagger}+b_{\alpha}\right)+\sum_{\alpha} \hbar \omega_{\alpha \sigma} b_{\alpha}^{\dagger} b_{\alpha}+H_{a n}\right)+t \sum_{\sigma}\right| \sigma\right\rangle\langle-\sigma|$

- All other index is other vibrational dof .  Which is not shifted. 
  
  For overlap between two state $|n^a\rangle =|n^{a}_{0} , n^{a}_1 , n^{a}_{2} , \cdots  \rangle $   , as $n_{0}^{a}$ is index for electronic state, $n_{1}^{a}$ is index for first vibrational state , then $n_{2}^{a} , n_{3}^{a} j, \cdots $ is quantum number for all other vibrational state.  

If we want overlap between two state : $\langle n^{a} | t (|\sigma\rangle \langle -\sigma|) |n^{b} \rangle \neq 0$ ,   ( this is coupling between electronic state as shown in Hamiltonian) , we must have :

$n_{0}^{a} \neq n_{0}^{b}$ , $n_{m}^{a} = n_{m}^{b}  (m=2,3,\cdots)$ , and there is extra factor as :

$\langle n_{1}^{a} | n_{1}^{b} \rangle $  : as first vibrational mode is shifted, we do not require $n_{1}^{a} = n_{1}^{b}$

#### Code for vib state $|m\rangle$ , $|n\rangle$ of first vib mode

Now about following code : 

```cpp
    # See Logan's note eq.(47).
     double Crossing_point_quanta_spin_down = pow(  (mfreq[1][0])/ (coupling_strength_to_mode0[0] + coupling_strength_to_mode0[1])
                                                   + 2 * coupling_strength_to_mode0[0] / (mfreq[0][1]) , 2 ) / 4 ;

    double Crossing_point_quanta_spin_up =
            pow(  (mfreq[1][0])/ (coupling_strength_to_mode0[0] + coupling_strength_to_mode0[1])
                  - 2 * coupling_strength_to_mode0[1] / (mfreq[0][1]) , 2 ) / 4  ;
```

eq.(47):

$m=\left[\frac{\epsilon}{4 \lambda}+\frac{\lambda}{\hbar \omega_{0}}\right]^{2}, \quad n=\left[\frac{\epsilon}{4 \lambda}-\frac{\lambda}{\hbar \omega_{0}}\right]^{2}$

**About frequency of first mode:**

mfreq[0][0] == mfreq[1][0] :  $\epsilon$  which is energy difference between electronic state.

```cpp
 double Crossing_point_quanta_spin_down = pow(  (mfreq[1][0])/ (coupling_strength_to_mode0[0] + coupling_strength_to_mode0[1])
                                                   + 2 * coupling_strength_to_mode0[0] / (mfreq[0][1]) , 2 ) / 4 ;
```

This equivalent to :

$ m = \frac{1}{4} [\frac{\epsilon}{\lambda_{1} + \lambda_{2} } + \frac{2 \lambda_{1}}{\hbar \omega_{0}}]^{2}$

Here :

 **coupling_strength_to_mode0[0] $\equiv \lambda_{1}$**    coupling strength of vib mode 0 to electronic state. for ground electronic state 

                              

  **coupling_strength_to_mode0[1]  $\equiv \lambda_{2}$**  oupling strength of vib mode 0 to electronic state. for excited electronic state 

```cpp
double Crossing_point_quanta_spin_up =
            pow(  (mfreq[1][0])/ (coupling_strength_to_mode0[0] + coupling_strength_to_mode0[1])
                  - 2 * coupling_strength_to_mode0[1] / (mfreq[0][1]) , 2 ) / 4  ;
```

This is equivalent to :

$n = \frac{1}{4} (\frac{\epsilon}{\lambda_{1} + \lambda_{2}} - \frac{2 \lambda_{2} }{\hbar \omega_{0} })^{2}$

( Here we assume $\lambda_{1}$ , $\lambda_{2}$ could be different  )

 **mfreq[1][0]** $\equiv \epsilon $  : vibrational ground state energy of excited electronic state

 **mfreq[0][1] = mfreq[1][1]** = $\omega_{0}$

```cpp
          crossing_region_state[0][0] = 0;
        crossing_region_state[0][1] = lround(Crossing_point_quanta_spin_down );

        crossing_region_state[1][0] = 1;
        crossing_region_state[1][1] = lround (Crossing_point_quanta_spin_up );Set initial state of system : 
```

- initial_vibrational_state[0]   : $|0 , m ,\cdots \rangle$ .  here 0 stands it is in ground electronic state.

$m$ stands its q.n. (quantum number) in first vibrational mode is == $m$  (see eq.(47) in Logan note and note above. )

- initial_vibrational_state[1] :  $|1  , n ,\cdots \rangle$ , here 1 stands for it is in excited electronic state. 

$n$ stands its q.n. in first vibrational mode is == n . (see eq.(47 ) in Logan's note).

#### 3. construct_Hamiltonian_basis_set( )

```cpp
void full_system:: construct_Hamiltonian_basis_set( ){
    // use this function we construct monomer state in a cube.
    if(my_id==0){
        int i, j, k;
        int i1;
        double  energy;

        double energy_for_vibrational_state;
        // ndetector0 and ndetector1 indicate current monomer mode quantum number (0001001 for example)  we are considering.
        vector<int> ndetector0(d.nmodes[0]);

        int location;
        bool exist=false;

        int state_one_norm_distance;

......
```

**Energy shift :**

```cpp
        // ground_electronic_state_initial_energy_shift : energy shift for ground state
        double ground_electronic_state_initial_energy_shift = 0 ;

        // excited_electronic_state_initial_energy_shift : energy shift for excited state.
        // this energy shift is due to different coupling strength in different electronic state. See eq.(22 b ) in Logan's note.
        double excited_electronic_state_initial_energy_shift = pow(coupling_strength_to_mode0[0],2) / d.mfreq[0][1] - pow(coupling_strength_to_mode0[1] , 2) / d.mfreq[0][1] ;
```

$\begin{aligned}
H_{0} &=\sum_{\sigma}|\sigma\rangle\langle\sigma|\left(\frac{\sigma}{2} \epsilon+\sigma \lambda_{\sigma}\left(b_{0}^{\dagger}+b_{0}\right)+\hbar \omega_{0 \sigma} b_{0}^{\dagger} b_{0}\right) \quad \equiv \sum_{\sigma} H_{0 \sigma} \\
&=\sum_{\sigma}|\sigma\rangle\langle\sigma|\left(\frac{\sigma}{2} \epsilon-\frac{\lambda_{\sigma}^{2}}{\hbar \omega_{0 \sigma}}+\hbar \omega_{0 \sigma} \tilde{b}_{0 \sigma}^{\dagger} \tilde{b}_{0 \sigma}\right)
\end{aligned}$

pow(coupling_strength_to_mode0[1] , 2) / d.mfreq[0][1]  : $\lambda_{1}^{2} / \omega_{0}$

pow(coupling_strength_to_mode0[0],2) / d.mfreq[0][1] : $\lambda_{0}^{2} / \omega_{0}$

Thus excited_electronic_state_initial_energy_shift = $-\lambda_{1}^{2}/\hbar \omega_{0} - (- \lambda_{0}^{2} / \hbar \omega_{0})$

**This is change of ground state energy due to shift** $b$ to $\tilde{b}$

**crossing region energy:**

```cpp
        double crossing_region_energy = d.mfreq[0][1]/2 *
                pow( d.mfreq[0][0] /(sqrt(2) *(coupling_strength_to_mode0[0] + coupling_strength_to_mode0[1]) ) +
                sqrt(2) * coupling_strength_to_mode0[0] / d.mfreq[0][1] , 2);
```

This is just $E = m \hbar \omega_{0} = n\hbar \omega_{0} + \epsilon$

#### Construct basis set for Hamiltonian:

```cpp
while (1) {
            label2:;  // label2 is for detector0 to jump to beginning of while(1) loop (this is inner layer of while(1))

            // ground electronic state
            ndetector0[0] = 0;
            // i1 begin from 1, this is because mode 0 is index to denote electronic state.
            // code below will let us go through all states (000000) - > (100000) -> (nmax0, nmax1, .... , nmax_n) with constraint (energy constrain + distance cutoff constraint) we add to sift states to use
            for (i1 = 1; i1 < d.nmodes[0]; i1++) {  // loop through detector0
                // define the way we loop through vibrational quantum numbe   (000000) - > (nmax0, nmax1, .... , nmax_n):
                ndetector0[i1] = ndetector0[i1] + 1;
                if (ndetector0[i1] <= d.nmax[0][i1]) break;
                if (ndetector0[d.nmodes[0] - 1] > d.nmax[0][d.nmodes[0] - 1]) {
                    ndetector0[d.nmodes[0] - 1] = 0;
                    goto label1;  // use goto to jump out of nested loop
                }
                ndetector0[i1] = 0;
.......
            }
```

This code  goes through all vibrational state at ground electronic state:

$|0 , \cdots \rangle $, covering all state whose quantum number < $n_{max}$ . We first increment index 1 , then index 2 , then index 3.

$|0,0,0, \cdots \rangle  \rightarrow |0,1,0,\cdots \rangle \rightarrow \cdots \rightarrow  |0, n_{max} , 0, \cdots \rangle \rightarrow |0,0,1,\cdots\rangle \rightarrow  \cdots \rightarrow |0,n_{max},1,\cdots\rangle $

```cpp
            energy = 0;
            energy_for_vibrational_state = 0;
            // calculate monomer 1 energy
            for (i = 0; i < d.nmodes[1]; i++) {
                energy = energy + ndetector0[i] * d.mfreq[1][i];
            }
            for( i = 1; i< d.nmodes[1]; i++){
                energy_for_vibrational_state = energy_for_vibrational_state + ndetector0[i] * d.mfreq[1][i];
            }
```

**energy** : this is energy for vibrational state + **ground state energy of electronic state**

**energy_for_vibrational_state** : this is energy for vibrational state.

```cpp
            //--------------------------------------------------------------------------------------------
            // criteria below make sure monomer 0 's energy is reasonable.
            // d.initial_vibrational_energy[0] only include vibrational energy of electronic state.
            if ( energy_for_vibrational_state > d.initial_vibrational_energy[1] + d.vibrational_energy_window_size) {
                // monomer 0's energy can not be larger than its initial energy + photon energy
                // jump to next monomer state.

                // start with first vibrational mode
                k = 1;
                while (ndetector0[k] == 0) {
                    ndetector0[k] = d.nmax[0][k];
                    k++;
                    if (k >= d.nmodes[0]) {
                        break;
                    }
                }
                if (k < d.nmodes[0]) {
                    ndetector0[k] = d.nmax[0][k];
                }
                goto label4;
            }
```

This will speed up process of going through vibrational state space. 

Because of the sequence we choose to go through state in state space, once energy is higher , then we can set all nonzero index to 0 and next 0 index to 1, for example :

$|1,2,3,0,1 , \cdots \rangle \rightarrow |0,0,0,1,1,\cdots \rangle  $ 

```cpp
            // criteria below make sure monomer 1 can not be too far away from bright state and lower bright state.
            // this distance only apply to vibrational state space.
            state_one_norm_distance = state_distance(ndetector0, d.initial_vibrational_state[1], d.nmodes[1]);
            if ( state_one_norm_distance > Rmax) {
                goto label4;
            }
```

We require quantum state we incorporate satisfy $| n - n_{init} | < R_{max} $

```cpp
location=find_position_for_insert_binary(monomer_qn_list0, ndetector0, exist);  // we check if this mode exist and the location we have to insert this state at the same time.
            if (!exist) {
                // when we push back we should consider arrange them in order. We compute location to insert in find_position_for_insert_binary() function:
                monomer_qn_list0.insert(monomer_qn_list0.begin() + location, ndetector0);
                dmat0.insert(dmat0.begin() + location, energy);
            }
```

To facilitate future reference of state, we also order state according to their quantum number. 

**find_position_for_insert_binary** is used to order the state according to q.n.

#### 

#### 4. construct_Hamiltonian_matrix()

```cpp
void monomer:: construct_Hamiltonian_matrix(ifstream & input, ofstream & output, ofstream & log,  vector<double> & dmat0,  vector<double> & dmat1,  vector<vector<int>> & monomer_qn_list0, vector<vector<int>> & monomer_qn_list1) {
    int m;

    // previously information for state space is in process 0. Now we broadcast this information to all procerss.
    construct_dv_dirow_dicol_dmatrix_MPI(log, dmat0, dmat1, monomer_qn_list0, monomer_qn_list1);

    // compute index for initial state (here state at crossing region) in Hamiltonian.
    compute_initial_vibrational_state_index();

    compute_monomer_offdiag_part_MPI(log,dmat0,dmat1,monomer_qn_list0,monomer_qn_list1);
 ....
```

##### 4.1 compute_initial_vibrational_state_index()

We know in some cases, the initial state we choose in program may be important for evaluating some quantity, thus we compute state index in state list according to q.n. of initial state.

##### 4.2 compute off diagonal matrix of Hamiltonian

```cpp
void monomer::compute_monomer_offdiag_part_MPI(ofstream & log,vector<double> & dmat0,  vector<double> & dmat1,  vector<vector<int>> & monomer_qn_list0, vector<vector<int>> & monomer_qn_list1)
```

**Coupling between state in different vibrational state**

**ground state for shifted harmonic oscillator** : 

$\begin{aligned}
|\alpha, 0\rangle=D(\alpha)|0\rangle &=e^{-\frac{1}{2} \alpha^{2}} e^{\alpha b^{\dagger}} e^{-\alpha b}|0\rangle \\
&=e^{-\frac{1}{2} \alpha^{2}} e^{\alpha b^{\dagger}}|0\rangle
\end{aligned}$

Here to compute coupling between vibrational state in different electronic state we have to consider overlap of two vibrational state : 

$t_{\left(\alpha_{\downarrow} m, \alpha_{\uparrow} n\right)}=t<\alpha_{\downarrow} m \mid \alpha_{\uparrow} n>$

$\begin{array}{r}
<\alpha_{\downarrow} m\left|\alpha_{\uparrow} n>=e^{-\frac{1}{2}\left(\alpha_{\downarrow}^{2}+\alpha_{\uparrow}^{2}\right)}<m\right| e^{-\alpha_{\downarrow} b^{\dagger}} e^{\alpha_{\downarrow} b} e^{\alpha_{\uparrow} b^{\dagger}} e^{-\alpha_{\uparrow} b} \mid n> \\
\left.=\exp \left(-\Delta \alpha^{2} / 2\right)<m\left|\exp \left(-\Delta \alpha b^{\dagger}\right) \exp (\Delta \alpha b)\right| n>\right)
\end{array}$

Here : $\Delta \alpha = \alpha _{\downarrow} - \alpha_{\uparrow}$

Above derivation use Baker–Campbell–Hausdorff formula :

[Baker–Campbell–Hausdorff formula - Wikiwand](https://www.wikiwand.com/en/Baker%E2%80%93Campbell%E2%80%93Hausdorff_formula)

More specifically : $e^{\alpha_{\downarrow} b} e^{\alpha_{\uparrow} b^{\dagger}} = e^{\alpha_{\downarrow}b + \alpha_{\uparrow}b^{\dagger}+\frac{1}{2} \alpha_{\downarrow} \alpha_{\uparrow}}$

$ e^{\alpha_{\uparrow} b^{\dagger}} e^{\alpha_{\downarrow} b}= e^{\alpha_{\downarrow}b + \alpha_{\uparrow}b^{\dagger}-\frac{1}{2} \alpha_{\downarrow} \alpha_{\uparrow}}$

```cpp
    double Alpha = (coupling_strength_to_mode0[0] + coupling_strength_to_mode0[1]) / (mfreq[0][1]);
```

Here $\alpha_{\downarrow} = \lambda_{1} / \hbar \omega_{0}$  , $\alpha_{\uparrow} = - \lambda_{2} / \hbar \omega_{0}$

thus $\Delta \alpha = (\lambda_{1} + \lambda_{2} ) / (\hbar \omega_{0} )$

We define $\Delta \alpha$ in code as    **Alpha** .

```cpp
           // Coupling within same electronic state.
            if((*vmode_ptr)[i][0] == (*vmode_ptr)[j][0] ){
                ntot=0;
                // k begins with 1, because index 0 is to denote different electronic state.
                for(k=1;k<nmodes[0];k++){
                    deln[k]= abs( (*vmode_ptr)[i][k] - (*vmode_ptr)[j][k] ); //  deln[k] = abs(monomer_vibrational_states_quantum_number_list[m][i][k] - monomer_vibrational_states_quantum_number_list[m][j][k]);
                    nbar[k]= sqrt(sqrt(double(max(1, (*vmode_ptr)[i][k])) * double(max(1, (*vmode_ptr)[j][k]  ))  )); // sqrt( vi[k] * vj[k] ).  This account for coefficient \sqrt{n + 1} generated by a^{+}|n> etc.
                    ntot= ntot+ deln[k];
                }
                // this is because we don't have q1 or q1 * q2 in Harmonic oscillator's Hamiltonian
                if (ntot == 2) {
                    for (k = 1; k < nmodes[0]; k++) {
                        if (deln[k] == 2) deln[k] = 4;
                        if (deln[k] == 1) deln[k] = 2;
                    }
                } else if (ntot == 1) {
                    for (k = 1; k < nmodes[0]; k++) {
                        if (deln[k] == 1) deln[k] = 3;
                    }
                } else if (ntot == 0) {
                    log << "Error! The off-diagonal element must differ in q.n." << endl;
                    MPI_Abort(MPI_COMM_WORLD,-8);
                }
                // check ntot with maxdis in quantum number space:
                // maxdis: 1-norm distance cutoff in state space for coupling between states.
                if (ntot < maxdis) {

                    if (ntot % 2 == 0) {
                        value = V_intra;  // V=0.03 as requirement.
                    } else {
                        value = -V_intra;
                    }

                    for (k = 1 ; k < nmodes[0]; k++) {
                        // aij is scaling factor for mode k.
                        value = value * pow(aij[0][k]* nbar[k], deln[k]);
                    }
                    // lij == V / (\Delta E) . Here V is anharmonic coupling between states, \Delta E is energy difference between states. Recall time-independent perturbation in QM.
                    // only when lij > cutoff, this anharmonic coupling will be important for dynamics.
                    if ( (*dmat_ptr)[i] != (*dmat_ptr)[j] ) {
                        lij = abs(value / ((*dmat_ptr)[i] - (*dmat_ptr)[j]));
                        if (lij > cutoff) {
                            monomer_mat[0].push_back(value);
                            monomer_irow[0].push_back(i);
                            monomer_icol[0].push_back(j);
                        }
                    } else {
                        monomer_mat[0].push_back(value);
                        monomer_irow[0].push_back(i);
                        monomer_icol[0].push_back(j);
                    }
                }
            }
```

Code for constructing off-diagonal matrix for coupling within same electronic state : 

This is just anharmonic vibrational coupling :

$V^{(\Delta n)}=(-1)^{\Delta n} V_{0} \times \prod_{i}\left(\left(\omega_{i} / \Omega\right)^{1 / 2} \times a\right)^{\Delta n_{i}}$

Or in Logan's note, that is eq.(15b , 15c):

$\begin{aligned}
&H_{0}=\sum_{\sigma}|\sigma\rangle\langle\sigma|\left(\frac{\sigma}{2} \epsilon-\frac{\lambda_{\sigma}^{2}}{\hbar \omega_{0 \sigma}}+\hbar \omega_{0 \sigma} \tilde{b}_{0 \sigma}^{\dagger} \tilde{b}_{0 \sigma}+\sum_{\alpha} \epsilon\left(\hat{n}_{\alpha}\right)\right) \quad \equiv \sum_{\sigma} H_{0 \sigma} \\
&H_{\phi}=\sum_{\sigma}|\sigma\rangle\langle\sigma|\left(\frac{1}{2} \sum_{\beta, \gamma}^{\prime} \phi_{0 \beta \gamma}\left(\tilde{b}_{0 \sigma}^{\dagger}+\tilde{b}_{0 \sigma}\right)\left(b_{\beta}^{\dagger}+b_{\beta}\right)\left(b_{\gamma}^{\dagger}+b_{\gamma}\right)+H_{\phi L W}\right) \quad \equiv \sum_{\sigma} H_{\phi \sigma}
\end{aligned}$

Code for coupling between electronic state : 

```cpp
            // Coupling between states with different electronic state.
            // Code below is for coupling between vibrational states with different electronic state.
            if( (*vmode_ptr)[i][0] != (*vmode_ptr)[j][0] ){
                // Check other vibrational state except vibrational mode with index: 1 is the same
                // as we mentioned before : all vib state except vib mode 0 (index1) should be equal to each other
                Other_vibrational_state_same = true ;
                for( k = 2 ; k <nmodes[0]; k ++  ) {
                    if ((*vmode_ptr)[i][k] != (*vmode_ptr)[j][k]) {
                        Other_vibrational_state_same = false;
                        break;
                    }
                }
                if( not Other_vibrational_state_same){
                    continue;
                }


                // include valid coupling
                // below we compute <spin_down, mode_1_index (m_index) | spin_up, mode_1_index (n_index) >
                // for <spin_up,  m_index , spin_down, n_index > == < spin_down, n_index | spin_up, m_index>
                if((*vmode_ptr)[i][0] == 0 and (*vmode_ptr)[j][0] == 1){
                    //  <spin_down | spin_up> case
                    m_index = (*vmode_ptr)[i][1];
                    n_index = (*vmode_ptr)[j][1];
                }
                else{
                    // <spin_up | spin_down> case:
                    m_index = (*vmode_ptr)[j][1];
                    n_index = (*vmode_ptr)[i][1];
                }
                // below compute <spin_down, m_index | spin_up, n_index >
                Minimum_Mode_1_index = min(  m_index , n_index ) ;
                tunneling_matrix_element = 0;

                for(k=0; k<= Minimum_Mode_1_index; k++){
                    // (-1)^{m_index - k}/ [(m_index - k)! * (n_index-k)! ] * sqrt(m! * n!) / (k!)
                    factorial_coefficient = 1 / (factorial( m_index - k ) * factorial(n_index - k) ) *
                            pow(-1 , m_index - k) * sqrt(factorial(m_index) * factorial(n_index) )/ (factorial(k));

                    tunneling_matrix_element = tunneling_matrix_element + pow(Alpha, m_index + n_index - 2 * k) *
                            factorial_coefficient;
                }
                // Coupling_between_electronic_state : t in Logan's note.
                tunneling_matrix_element = tunneling_matrix_element * exp(- pow(Alpha,2) / 2 ) * Coupling_between_electronic_state;

                if( i == initial_state_index[0] and j== initial_state_index[1]){
                    cout << "Found " << endl;
                }

                // cutoff criteria. This equivalent to strong coupling is among states in two electronic state which have similar energy. (near crossing region)


                if ( (*dmat_ptr)[i] != (*dmat_ptr)[j] ) {
                    lij = abs(tunneling_matrix_element / ((*dmat_ptr)[i] - (*dmat_ptr)[j]));
                    if (lij > cutoff_for_coupling_between_electronic_state) {
                        monomer_mat[0].push_back( tunneling_matrix_element );
                        monomer_irow[0].push_back(i);
                        monomer_icol[0].push_back(j);
                    }
                } else {
                    monomer_mat[0].push_back( tunneling_matrix_element );
                    monomer_irow[0].push_back(i);
                    monomer_icol[0].push_back(j);
                }

            }
```

To compute : 

$\begin{array}{r}
<\alpha_{\downarrow} m\left|\alpha_{\uparrow} n>=e^{-\frac{1}{2}\left(\alpha_{\downarrow}^{2}+\alpha_{\uparrow}^{2}\right)}<m\right| e^{-\alpha_{\downarrow} b^{\dagger}} e^{\alpha_{\downarrow} b} e^{\alpha_{\uparrow} b^{\dagger}} e^{-\alpha_{\uparrow} b} \mid n> \\
\left.=\exp \left(-\Delta \alpha^{2} / 2\right)<m\left|\exp \left(-\Delta \alpha b^{\dagger}\right) \exp (\Delta \alpha b)\right| n>\right)
\end{array}$

We use following code : 

```cpp
                // below compute <spin_down, m_index | spin_up, n_index >
                Minimum_Mode_1_index = min(  m_index , n_index ) ;
                tunneling_matrix_element = 0;

                for(k=0; k<= Minimum_Mode_1_index; k++){
                    // (-1)^{m_index - k}/ [(m_index - k)! * (n_index-k)! ] * sqrt(m! * n!) / (k!)
                    factorial_coefficient = 1 / (factorial( m_index - k ) * factorial(n_index - k) ) *
                            pow(-1 , m_index - k) * sqrt(factorial(m_index) * factorial(n_index) )/ (factorial(k));

                    tunneling_matrix_element = tunneling_matrix_element + pow(Alpha, m_index + n_index - 2 * k) *
                            factorial_coefficient;
                }
                // Coupling_between_electronic_state : t in Logan's note.
                tunneling_matrix_element = tunneling_matrix_element * exp(- pow(Alpha,2) / 2 ) * Coupling_between_electronic_state;
```

**Explanation :** 

$\langle m | \exp(-\Delta \alpha \hat{b}^{\dagger}) \exp(\Delta \alpha \hat{b} )|n\rangle  = \\ \langle m| (1 - \Delta \alpha \hat{b}^{\dagger} + \frac{(-\Delta\alpha)^{2} (\hat{b}^{\dagger})}{2!} + \cdots )  (1 + \Delta \alpha \hat{b} + \cdots ) |n\rangle  \\\ = \sum_{k} \frac{(- \Delta \alpha)^{m-k} }{(m-k)!} [\langle m| (\hat{b}^{\dagger})^{m-k}]  \times \frac{ (\Delta \alpha)^{n-k} }{(n-k)!} (\hat{b}^{n-k} |n \rangle ) \\ = \sum_{k} \frac{(- \Delta \alpha)^{m-k} }{(m-k)!} [\langle k| \frac{\sqrt{m!}}{\sqrt{k!}}  ] \times \frac{ (\Delta \alpha)^{n-k} }{(n-k)!} ( \sqrt{\frac{n!}{k!}} |k \rangle ) =\\ \sum_{k=0}^{min(m,n)} (\Delta \alpha)^{m+n-2k} \frac{(-1)^{m-k}}{(m-k)! (n-k)!} \frac{\sqrt{m! n!}}{k!}$ 

We see this is exactly what code above compute.

##### 4.3 construct sampling state index

```cpp
    if( not load_sampling_state ){
        // construct sampling_state_index which contain states we want to compute IPR
        // we divide all state to 20 categories according to energy, and sample 10 states in each cateogroty
        construct_sampling_state_index(dmat0);
    }
    else{
        // load sampling_state_index from file: sampling_state_info
        load_sampling_state_index(log);
    }
```

We choose state to do simulation.

### Evolve_full_sysem_MPI()

Evolve Schrodinger equation

$i\hbar \frac{d \psi}{dt} = H \psi $

#### 

#### 1. d.prepare_evolution():

<img title="" src="file:///home/phyzch/CLionProjects/LW_code_clean_version/note%20fig/PNG%20image.png" alt="" width="596">

![](/home/phyzch/CLionProjects/LW_code_clean_version/note%20fig/PNG%20image1.png)

We use MPI to speed up matrix array multiplication. 

We divide wave function $\psi$ into different part and store in different processes.

$d\psi = \psi(t+dt) - \psi = dt \times (H \psi)$

We know Hamiltonian $H$ as a matrix may have component : $H(m,n)\neq 0$ , and element $m,n$ belong to different process

Hamiltonian matrix $H$ is shared in different processes.

For process $p_n$, it has array of $\psi$ from index $n_{1} \sim n_{N}$ , if we have $H(n_{j},m)$ where $n_{j} \in [n_{1} , n_{N}]$ but $m \not \in [n_{1} , n_{N}]$ .

 Say we have in total $M$ element that $H(n_{j},m)\neq 0$ , then we  prolong the array that store $\psi$ from size $n_{N} - n_{1}$ to $n_{N} - n_{1} + M$ . with local array of extra size $M$.

Then each time before we compute $H \psi$ , we send data $\psi[m]$ from other processes and store in extra M elements in local array. Then when compute $H(n_{j},m)$  we multiply H by element in local array.

```cpp
    to_recv_buffer_len[0] = construct_receive_buffer_index(remoteVecCount[0],remoteVecPtr[0],
            remoteVecIndex[0],0);  // construct buffer to receive.
    to_send_buffer_len[0]= construct_send_buffer_index(remoteVecCount[0],remoteVecPtr[0],remoteVecIndex[0],
                                                    tosendVecCount[0],tosendVecPtr[0], tosendVecIndex[0]);
```

This code analyze Hamiltonian H and find all the index for $m$ of  vector  $\psi$ that should be trasferred between different process before matrix -vector multiplication.

#### 1.1 construct_receive_buffer_index()

```cpp
vector <int> col_index_copy = monomer_icol[detector_index];  // column index in Hamiltonian has nonzero element.
    sort(col_index_copy.begin(),col_index_copy.end()); // sort vector.
    int col_array_size = col_index_copy.size();
    int prev_col=-1;
    j=0;
    for(i=0;i<col_array_size;i++){
        if( (col_index_copy[i]>prev_col)   and ( (col_index_copy[i]<local_begin) or (col_index_copy[i] >= local_end) )  ){
            // this matrix element is not in process.
            if (col_index_copy[i] >= vsize * (num_proc-1) ){
                remote_pc_id = num_proc-1;
            }
            else{
                remote_pc_id = col_index_copy[i] / vsize;
            }
            remoteVecCount_element[remote_pc_id] ++;
            remoteVecIndex_element [j] = col_index_copy[i];  // vector index need to receive. (global index , ordered)
            j++;
        }
        prev_col= col_index_copy[i];
    }
    remoteVecPtr_element[0]=0;   // displacement for remote vector from each process in remoteVecIndex.
    for(i=1;i<num_proc;i++){
        remoteVecPtr_element[i] = remoteVecPtr_element[i-1] + remoteVecCount_element[i-1];
    }
    for(i=0;i<num_proc;i++){
        total_remoteVecCount = total_remoteVecCount + remoteVecCount_element[i];
    }
```

Here **monomer_icol[detector_index]** is column index for nonzero element in Hamiltonian reside in given process.

1 . We record index in  remoteVecIndex_element

2. We record number of array index need to transfer in remoteVecCount_element

3. remoteVecPtr_element : used in MPI. Tell MPI the displacement of vectors.
   
   We use MPI_all_to_allv

#### 1.2 constuct_send_buffer_index()

Prepare index information to send buffer to other process.

use MPI_Alltoall and MPI_Alltoallv 

#### 1.3 Allocate space for element received from other process.

```cpp
    for(m=0;m< sampling_state_list_size;m++){
        xd[m].resize(monomer_matsize[0] + to_recv_buffer_len[0]);
        yd[m].resize(monomer_matsize[0] + to_recv_buffer_len[0]);
        recv_xd[m] = new double [to_recv_buffer_len[0]];
        recv_yd[m]= new double [to_recv_buffer_len[0]];
        send_xd[m]= new double [to_send_buffer_len[0]];
        send_yd[m] = new double [to_send_buffer_len[0]];
    }
```

This code extend $xd$ (real part of wave function) , $yd$ (imaginary part of wave function)  array to reside extra M elements. 

#### 1.4 Decide row, index array for matrix multiplication

```cpp
    // construct local_dirow, local_dicol
    local_dirow[0].reserve(monomer_matnum[0]);
    local_dicol[0].reserve(monomer_matnum[0]);
    for(i=0;i<monomer_matnum[0];i++){
        local_dirow[0].push_back(monomer_irow[0][i] - my_id * vsize);  // set local index for row index
        col_index_to_search= monomer_icol[0][i];
        search_Ind=(int *) bsearch(&col_index_to_search,remoteVecIndex[0],to_recv_buffer_len[0],sizeof(int),compar);
        if(search_Ind!=NULL){
            // this column index is not in local matrix, and we should get it from other process (remoteVec)
            local_dicol[0].push_back(monomer_matsize[0] + (search_Ind-remoteVecIndex[0]) );
        }
        else{ // this column index is in local matrix.
            local_dicol[0].push_back (monomer_icol[0][i] - my_id * vsize );
        }
    }
```

As in Matrix Multiplication, we may use array element received from other process, we have to construct separate irow, icol for matrix vector multiplication.

Here remoteVecIndex is index in rang $[n_{1} , n_{N}]$ for process $p_{n}$ . 

local_dirow , local_dicol is used in matrix vector multiplication as row and column index for Hamiltonian H.

$d\psi [row] = H [row , col] \times \psi[col]$

local_dirow is always element in range $[n_{1} , n_{N}]$ for process $p_{n}$ ,

local_dicol  is in $[n_{1} , n_{N}]$ if $H[n,m]$ has $m \in [n_{1} , n_{N}]$ , otherwise local_dicol is in range $[n_{N} , n_{N} + M]$ .

#### 2 . Evolve wave function $\psi$       SUR_one_step():

```cpp
void monomer::SUR_onestep_MPI(){
    int m, i;
    int irow,icol;
    int nearby_state_list_size = sampling_state_index.size();
    // do simulation starting from different state

    // update imaginary part of wave function calculated in different process.
    update_dy(nearby_state_list_size);
    // SUR algorithm
    for(i=0;i<monomer_matnum[0];i++){
        // make sure when we compute off-diagonal matrix, we record both symmetric and asymmetric part
        irow = local_dirow[0][i];
        icol = local_dicol[0][i]; // compute to point to colindex in
        for(m=0;m<nearby_state_list_size;m++) {
            xd[m][irow] = xd[m][irow] + monomer_mat[0][i] * yd[m][icol] * cf;
        }
    }

    // update real part of wave function calculated in different process.
    update_dx(nearby_state_list_size);
    // SUR algorithm.
    for(i=0;i<monomer_matnum[0];i++){
        irow= local_dirow[0][i];
        icol= local_dicol[0][i];
        for(m=0;m<nearby_state_list_size;m++){
            yd[m][irow] = yd[m][irow] - monomer_mat[0][i] * xd[m][icol] * cf;
        }
    }

}
```

See before we do  $H \times \psi$ we have  : **update_dx()**  , **update_dy()** .

These functions are for upodating elements of $\psi(m)$ for $H(n,m) \neq 0$  if $m$ is not in local array of wave function $\psi$ .
