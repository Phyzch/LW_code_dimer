//
// Created by phyzch on 4/15/20.
// This part of code estimate the cost of resources in our program
//
#include"system.h"
#include"util.h"
using namespace std;
typedef struct{
    unsigned long size, resident,share,text,lib,data,dt;
}statm_t;

