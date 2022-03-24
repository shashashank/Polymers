#include<cstdlib>
#include<cstdio>
#include<cmath>
#include<iostream>
#include<fstream>
#include<vector>
#include<random>
#include<ctime>
#include<dirent.h>
#include<string> //for string
#include<sstream> //covert into to string
#include <sys/types.h> //mkdir
#include <sys/stat.h>   //mkdir
#include <iomanip> // setprecision
#include <experimental/filesystem>
#include <omp.h>

using namespace std;

int readFile(string file, vector<string>& data);
void extractConfig(string oldFile, vector<double>& X, vector<double>& Y, int R);

