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

int readFile(std::string file, std::vector<std::string>& data);
void extractConfig(std::string oldFile, std::vector<double>& X, std::vector<double>& Y, int R);

