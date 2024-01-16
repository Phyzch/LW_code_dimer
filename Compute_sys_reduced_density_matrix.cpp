//
// Created by phyzch on 6/20/20.
// This file include function to compute system reduced density matrix fast.
//

#include<iostream>
#include"system.h"
#include"util.h"
using namespace std;

int compare_vmode (const vector<int> & vmode1, const  vector<int> & vmode2, const vector<int> & vmode1_0, const vector<int> & vmode2_0){
    // first compare monomer_qn_list1, in this way, we don't have to merge sort between different process.
    int size1= vmode1.size();
    if( size1 != vmode1_0.size()){
        cout<< "error! two vector's mode should equal."<<endl;
        exit(-8);
    }
    int size2= vmode2.size();
    if(size2 != vmode2_0.size()){
        cout<< "error! two vector's mode should equal."<<endl;
        exit(-8);
    }

    int i;
    for(i=0;i<size1;i++){
        if(vmode1[i]>vmode1_0[i]) return 1;
        else if (vmode1[i] < vmode1_0[i]) return -1;
    }
    for(i=0;i<size2;i++){
        if(vmode2[i]>vmode2_0[i]) return 1;
        else if (vmode2[i]<vmode2_0[i]) return -1;
     }
    return 0;
}


// ---------- end of merge_sort algorithm ----------------------------------------------
