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

#include "poly_param.c"
#include "extra.h"

constexpr double pi = 22.000 / 7;

mt19937 generator(time(NULL));
normal_distribution<double> gau_dist(0.0, 1.0);

vector<double> Px(N), Py(N);
vector<double> Fx(N), Fy(N);
vector<double> intrX(N), intrY(N);
vector<omp_lock_t> iXLock(N), iYLock(N);
double dt = 1e-3;
double radius = 1.5;
double l = 1.0;
int k = 100;
double sqrt_dt = sqrt(dt);
double MAXIT =  par_maxIT / dt;
double start_time;


double spring(double x, double y, double length, int spring_constant);
void animate(void);
void initialize(void);
void write_VMD_data(ostream& os, int it);
void noise(void);
double e2e_distance(double x_ini, double x_fin, double y_ini, double y_fin);
void SForce(int P1, int P2);
void IForce(int P1, int P2);

int main(int argc, char *argv[])
{
	start_time = omp_get_wtime();
	ofstream out("data_poly.xyz");
	ofstream o("vmd_data_poly.xyz");
	string old = "old.xyz";

	int skip = 0;
	if (std::experimental::filesystem::v1::exists(old)){
		skip = 1;
	  	extractConfig(old, Px, Py);
	}
	else initialize();

	out << "  no of particles::" << N << "  length: l:: " << l << endl;

	double frame = par_maxFrame;
	int tn = MAXIT / (frame * 10);
	double e2eVal = e2e_distance(Px[0], Px[N-1], Py[0], Py[N-1]);

	if (skip==0){
		for (int it = 0; it < MAXIT; it++){
			animate();
		}
	}

#pragma omp parallel for
	for (int ilock = 0; ilock<N; ilock++){
		omp_init_lock(&iXLock[ilock]);
		omp_init_lock(&iYLock[ilock]);
	}

	for (int it = 0; it < MAXIT; it++){
		if ((it % tn == 0) && (writeFlag==1)){
			write_VMD_data(o, it);
		}
		animate();
	}

#pragma omp parallel for
	for (int ilock = 0; ilock<N; ilock++){
		omp_destroy_lock(&iXLock[ilock]);
		omp_destroy_lock(&iYLock[ilock]);
	}

	out << "Time elapsed: " << ( (omp_get_wtime() - start_time) / (double) CLOCKS_PER_SEC) << endl;
}


void initialize(void){
	for (int i=0; i<N; i++)
	{
		Px[i] = double(i * l);
		Py[i] = 0.0;
	}
}

void animate(void){
	double Fp, Fn, Fxi, Fyi, Fsi, Fsj;
	double xp, xn, xsi, xsj, yp, yn, ysi, ysj;
	double Pxi, Pyi;
	double phip, phin, phisi, phisj;

	omp_set_dynamic(1);
	omp_set_num_threads(omp_get_max_threads());

// add schedule(guided, 10) for longer chains
#pragma omp parallel for
	for (int i=0; i<N; i++)
	{
	  	Fx[i] = 0.0; Fy[i] = 0.0;
	  	if (i!=0){
			SForce(i, i-1); //spring force from previous particle in chain
	  	}
	  	if (i!=(N-1)){
			SForce(i, i+1); //spring force from next particle in chain
		}
		for (int j=i+2; j<N; j++){
			double apart = e2e_distance(Px[i], Px[j], Py[i], Py[j])+ intrFlag;
			if (apart<radius){
				IForce(i, j);
			}
		}
	}
// add schedule(static, 100) for longer chains
#pragma omp parallel for 
	for (int i=0; i<N; i++){
		Px[i] += (Fx[i] + intrX[i])*dt + gau_dist(generator)*sqrt_dt;
		Py[i] += (Fy[i] + intrY[i])*dt + gau_dist(generator)*sqrt_dt;
		intrX[i] = 0.0;
		intrY[i] = 0.0;
	}
}


double spring(double x, double y, double length, int spring_constant) {
    k = double(spring_constant);
    double distance = sqrt(x*x + y*y);
    return k*(distance - length);
}

//Spring Force is calcuated and handed to the force vectors
void SForce(int P1, int P2){
	double x, y, phi, F;
	x = Px[P2] - Px[P1];
	y = Py[P2] - Py[P1];
	phi = atan2(y, x);
	F = spring(x, y, l, k);
	Fx[P1] += F*cos(phi);
	Fy[P1] += F*sin(phi);
}


//Interaction Force
void IForce(int P1, int P2){
	double x, y, phi1, phi2, F;
	x = Px[P2] - Px[P1];
	y = Py[P2] - Py[P1];
	F = spring(x, y, radius, k);
	phi1 = atan2(y, x);
	phi2 = atan2((-1.0*y), (-1.0*x));

	omp_set_lock(&iXLock[P1]);
	intrX[P1] += F*cos(phi1);
	omp_unset_lock(&iXLock[P1]);

	omp_set_lock(&iYLock[P1]);
	intrY[P1] += F*sin(phi1);
	omp_unset_lock(&iYLock[P1]);

	omp_set_lock(&iXLock[P2]);
	intrX[P2] += F*cos(phi2);
	omp_unset_lock(&iXLock[P2]);

	omp_set_lock(&iYLock[P2]);
	intrY[P2] += F*sin(phi2);
	omp_unset_lock(&iYLock[P2]);
}

void write_VMD_data(ostream& os, int it)
{
	os << N << endl;
	os << it << endl;
	os.setf(ios::fixed, ios::floatfield);
	os.precision(5);
	for (int i = 0; i < N; i++)
	{
		os << "s " << Px[i] << " " << Py[i] << endl;
	}
	os.flush();
}

double e2e_distance(double x_ini, double x_fin, double y_ini, double y_fin)
{
	double X = x_fin - x_ini;
	double Y = y_fin - y_ini;
	return sqrt(X*X + Y*Y);
}