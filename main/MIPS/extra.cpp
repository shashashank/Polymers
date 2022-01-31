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

void extractConfig(string oldFile, string newFile, int dim){
  ifstream infile(oldFile);
  ofstream outfile(newFile);
  vector<string> fileLines;
  string line;
  while (getline(infile, line)){
    fileLines.push_back(line);
  }
  istringstream iss(fileLines[0]);
  int N;
  iss >> N;
  printf("Total number of particles: %i \n", N);
  int fileLength = size(fileLines);
  int iters = fileLength / (N + 2);
  int t_pos = (iters - 1) * (N + 2);
  char s;
  double x, y, z;
  int count = 0;
  if (dim==2){
    for (int i = t_pos + 2; i < fileLength; ++i){
      istringstream tmp(fileLines[i]);
      tmp >> s >> x >> y;
      outfile << x << "\t" << y << endl;
      ++count;
    }
  }
  else if (dim==3){
    for (int i = t_pos + 2; i < fileLength; ++i){
      istringstream tmp(fileLines[i]);
      tmp >> s >> x >> y >> z;
      outfile << x << "\t" << y << "\t" << z << endl;
      ++count;
    }
  }
  printf("Number of particles transcribed: %i \n", count);
}

void importConfig(string config, vector<double>& X, vector<double>& Y, vector<double>& THETA){
  ifstream infile(config);
  vector<string> fileLines;
  string line;
  int i = 0;
  double a, b, c;
  while (getline(infile, line)){
    istringstream iss(line);
    iss >> a >> b >> c;
    X[i] = a;
    Y[i] = b;
    THETA[i] = c;
    ++i;
  }
}

// int main(){
//   extractConfig("vmd_data.xyz", "config", 3);
// }