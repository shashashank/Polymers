#include<cstdlib>
#include<cstdio>
#include<cmath>
#include<iostream>
#include<fstream>
#include<vector>
#include<random>
#include"strtk.hpp"
#include<ctime>
#include<dirent.h>
#include<string> //for string
#include<sstream> //covert into to string
#include <sys/types.h> //mkdir
#include <sys/stat.h>   //mkdir
#include <iomanip> // setprecision

#include "extra.h"
using namespace std;

void extractConfig(void);
void importConfig(vector<double>, vector<double>);

void extractConfig(){
  ifstream infile("vmd_data_poly_old.xyz");
  ofstream outfile("config");
  vector<string> fileLines;
  string line;
  while (getline(infile, line)){
    fileLines.push_back(line);
  }
  istringstream iss(fileLines[0]);
  int N;
  iss >> N;
  cout << N << endl;
  int fileLength = size(fileLines);
  int iters = fileLength / (N + 2);
  int t_pos = (iters - 1) * (N + 2);
  char s;
  double x, y;
  int count = 0;
  for (int i = t_pos + 2; i < fileLength; ++i){
    istringstream tmp(fileLines[i]);
    tmp >> s >> x >> y;
    outfile << x << "\t" << y << endl;
    ++count;
  }
  cout << count <<endl;
}

void importConfig(vector<double>& Px, vector<double>& Py){
  ifstream infile("config");
  vector<string> fileLines;
  string line;
  int i = 0;
  double x, y;
  while (getline(infile, line)){
    istringstream iss(line);
    iss >> x >> y;
    Px[i] = x;
    Py[i] = y;
    ++i;
  }
}

// std::string line;
// while (std::getline(infile, line))
// {
//     std::istringstream iss(line);
//     int a, b;
//     if (!(iss >> a >> b)) { break; } // error
//     // process pair (a,b)
// }
