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
#include <omp.h>

using namespace std;



vector<double> x(N), y(N);



int main(){
    initialize_square_lattice();
}

void initialize_square_lattice(){
    int nptl = 0;
	for (int i = 0; i < n; i++){
		for (int j = 0; j < n; j++){
			if (j % 2 == 0)
				x[nptl] = -L / 2.0 + i * a;
			else
				x[nptl] = -L / 2.0 + (2.0 * i + 1.0) * a / 2.0;
			if (i == 0)
				y[nptl] = -L / 2.0 + (2.0 * j + 1.0) * a / 2.0;
			else
				y[nptl] = -L / 2.0 + j * a;

			theta[nptl] = sqrt(2 * Dr * dt) * gau_dist(generator);
			nptl++;
		}
    }
}

