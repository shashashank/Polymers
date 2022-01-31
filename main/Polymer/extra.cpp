#include<cstdlib>
#include<cstdio>
#include<cmath>
#include<iostream>
#include<fstream>
#include<vector>
#include<random>
//#include"strtk.hpp"
#include<ctime>
#include<dirent.h>
#include<string> //for string
#include<sstream> //covert into to string
#include <sys/types.h> //mkdir
#include <sys/stat.h>   //mkdir
#include <iomanip> // setprecision

using namespace std;
#include "extra.h"

void extractConfig(string oldFile, vector<double>& X, vector<double>& Y, int R){
    ifstream infile(oldFile);
    vector<string> fileLines;
    string line;
    int N, t_pos, fileLength, iters;
    while (getline(infile, line)){
        fileLines.push_back(line);
    }
    istringstream iss(fileLines[0]);
    iss >> N;
    printf("Total number of particles: %i \n", N);
    fileLength = fileLines.size();
    iters = fileLength / (N + 2);
    if (R==0){
        t_pos = (iters - 1) * (N + 2);
    }else{
        random_device rd;
        mt19937 gen(rd());
        uniform_int_distribution<> distr(1, iters);
        t_pos = distr(gen) * (N + 2);
    }
    char s;
    double a, b;
    int count = 0;
    for (int i = t_pos + 2; i < (t_pos + 2 + N); ++i){
        istringstream iss(fileLines[i]);
        iss >> s >> a >> b;
        X[count] = a;
        Y[count] = b;
        ++count;
    }
    printf("Number of particles initialised: %i \n", count);
}
