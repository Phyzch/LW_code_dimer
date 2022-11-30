//
// Created by phyzch on 6/18/20.
//
#pragma once
#ifndef QUANTUM_MEASUREMENT_UTIL_H
#define QUANTUM_MEASUREMENT_UTIL_H

#include<iostream>
#include<time.h>
#include<stdio.h>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <sys/stat.h>
#include <experimental/filesystem>
#include<iomanip>
#include <complex>
#include <assert.h>
#include <vector>
#include<ctime>
#include<algorithm>
#include<stdlib.h>
#include<mpi/mpi.h>
#include<sys/resource.h>
//using namespace concurrency;
#define pi2 3.141592653589793*2
using namespace std;
extern double energy_window_size;
extern int Rmax;

extern bool Turn_on_Vanvleck;
extern int ndegre;
extern int ndegrx2;
extern double es_criteria;

extern double detector_coupling_time;
extern bool Turn_on_Gaussian_coupling;

// define function here
float ran2(long& idum);
void estimate_memory_cost(ofstream & resource_output);  // output resource cost to file at give time step.
void convert_dv(const vector<vector<int>> & vec_2d, vector <int>  & vec_1d , vector <int> & displacement , vector <int> & element_size );
// used for cnostruct buffer for communication between process for matrix multiplication.
int construct_send_buffer_index(int * remoteVecCount, int * remoteVecPtr, int * remoteVecIndex, int * tosendVecCount_element, int * tosendVecPtr_element, int * & tosendVecIndex_ptr);

int compar(const void * a, const void * b);

#endif //QUANTUM_MEASUREMENT_UTIL_H

