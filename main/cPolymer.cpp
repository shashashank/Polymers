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

clock_t start = clock();

constexpr double pi = 22.000 / 7;

mt19937 generator(time(NULL));
normal_distribution<double> gau_dist(0.0, 1.0);

vector<double> Px(N), Py(N);
vector<double> Fx(N), Fy(N);
vector<double> dDist(par_maxIT);
double dt = 1e-3;
double radius = 1.5;
double l = 1.0;
int k = 100;
double sqrt_dt = sqrt(dt);
double MAXIT =  par_maxIT / dt;


double spring(double x, double y, double length, int spring_constant);
double animate(void);
void initialize(void);
void write_VMD_data(ostream& os, int it);
void noise(void);
double e2e_distance(double x_ini, double x_fin, double y_ini, double y_fin);


int main(int argc, char *argv[])
{
	ofstream out("data_poly.xyz");
	ofstream ete("end-to-end.d");
	ofstream o("vmd_data_poly.xyz");

	int skip = 0;
	if (std::experimental::filesystem::v1::exists("vmd_data_poly_old.xyz")){
	  skip = 1;
	  extractConfig();
	  importConfig(Px, Py);
	}
	else initialize();

	out << "  no of particles::" << N << "  length: l:: " << l << endl;

	ete << "T" << "\t" << "D" << endl;

	double frame = par_maxFrame;
	int tn = MAXIT / (frame * 10);
	double e2eVal = e2e_distance(Px[0], Px[N-1], Py[0], Py[N-1]);

	if (skip==0){
	  for (int it = 0; it < MAXIT; it++){
	    animate();
	  }
	} else {
	  for (int it = 0; it < MAXIT/10; it++){
	    animate();
	  }
	}
        for (int it = 0; it < MAXIT/10; it++)
	{
	  if ((it % tn == 0) && (writeFlag==1)){
			write_VMD_data(o, it);
		}
       	       	ete << it + 1 << "\t" <<  e2eVal << endl;
		e2eVal = animate();
	}
	out << "Time elapsed: " << ( (double(clock()) - start) / CLOCKS_PER_SEC) << endl;
}


void initialize(void){
	for (int i=0; i<N; i++)
	{
		Px[i] = double(i * l);
		Py[i] = 0.0;
	}
	for (int i=0; i<par_maxIT; i++)
	{
		dDist[i] = 0.0;
	}
}

double animate(void){
	double Fp, Fn, Fxi, Fyi, Fsi, Fsj;
	double xp, xn, xsi, xsj, yp, yn, ysi, ysj;
	double Pxi, Pyi;
	double phip, phin, phisi, phisj;
	vector<double> intrX(N), intrY(N);

	Pxi = 0.0;
	Pyi = 0.0;
    omp_set_num_threads(8);
#pragma omp parallel for schedule(guided, 10)
	for (int i=0; i<N; i++)
	{
	  Fx[i] = 0.0; Fy[i] = 0.0;
	  if (i!=0){
	    xp = Px[i-1] - Px[i];
	    yp = Py[i-1] - Py[i];
	    phip = atan2(yp, xp);
	    Fp = spring(xp, yp, l, k);
	    Fx[i] -= Fp*cos(phip);
	    Fy[i] -= Fp*sin(phip);
	  }
	  if (i!=(N-1)){
	    xn = Px[i+1] - Px[i];
	    yn = Py[i+1] - Py[i];
	    phin = atan2(yn, xn);
	    Fn = spring(xn, yn, l, k);
	    Fx[i] -= Fn*cos(phin);
	    Fy[i] -= Fn*sin(phin);
	  }
	  for (int j=i+2; j<N; j++){
	    double apart = e2e_distance(Px[i], Px[j], Py[i], Py[j])+ intrFlag;
	    if (apart<radius){
	      xsi = Px[j] - Px[i];
	      ysi = Py[j] - Py[i];
	      Fsi = -spring(xsi, ysi, radius, k);
	      phisi = atan2(ysi, xsi);
	      intrX[i] += Fsi*cos(phisi);
	      intrY[i] += Fsi*sin(phisi);
	      xsj = Px[i] - Px[j];
	      ysj = Py[i] - Py[j];
	      Fsj = -spring(xsj, ysj, radius, k);
	      phisj = atan2(ysj, xsj);
	      intrX[j] += Fsj*cos(phisj);
	      intrY[j] += Fsj*sin(phisj);
	    }
	  }
	}
#pragma omp parallel for schedule(static, 100)
	for (int i=0; i<N; i++){
	  Px[i] += (Fx[i] + intrX[i])*dt + gau_dist(generator)*sqrt_dt;
	  Py[i] += (Fy[i] + intrY[i])*dt + gau_dist(generator)*sqrt_dt;
	  intrX[i] = 0.0;
	  intrY[i] = 0.0;
	}
	return e2e_distance(Px[0], Px[N-1], Py[0], Py[N-1]);
}


double spring(double x, double y, double length, int spring_constant) {
    k = double(spring_constant);
    double distance = sqrt(x*x + y*y);
    return -k*(distance - length);
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
