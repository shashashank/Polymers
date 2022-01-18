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
#include "param.c"
clock_t start = clock();

using namespace std;

int N1 = par_N;
double phi = par_phi;
double Pe = par_Pe;
double alpha = 6;
constexpr int sigma = 1;

constexpr double D = 1;
double Dr = 3 * D / (sigma*sigma);

constexpr double pi = 22.000 / 7;
/* time limit */
constexpr double TAU = 1;
constexpr double dt = 0.00001 * TAU;
double MAXIT =  par_maxIT * TAU / dt;
/*************/


/******calculate parameter for cell list*******/
double rcut = sigma * pow(2, 1.0 / alpha);
int n = sqrt(N1);
int N = n * n;
double L = sqrt(N * pi * sigma * sigma * 0.25 / phi);
int M = int(L / rcut);
double a = 2 * L / (2 * n + 3);

double BOXL = M * rcut;
int ncell = M * M;
int mapsiz = 4 * ncell;
double A = BOXL * BOXL;

/***********  ***********/

vector<int> list_(N);
vector<double> Fx(N), Fy(N);
vector<int> head(ncell + 1);
vector<int> map_(mapsiz + 2);
int ix, iy, iz, imap;

vector<double> theta(N);
vector<double> x(N), y(N);
vector<double> vx(N);
vector<double> vy(N);

mt19937 generator(time(NULL));
normal_distribution<double> gau_dist(0.0, 1.0);
uniform_real_distribution<double> uni_dist(-pi / 4, pi / 4);


/*******************sub function*****************/
void initialize_square_lattice()
{
	int nptl = 0;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
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

/***********cell list*************/
int icell(int ix, int iy)
{
	int a = 1 + ( ix - 1 + M) % M + (( iy - 1 + M) % M ) * M ;
	return a;
}
void maps()
{
	for (int iy = 1; iy <= M; iy++)
	{
		for (int ix = 1; ix <= M; ix++)
		{
			imap = ( icell ( ix, iy) - 1 ) * 4;
			map_[ imap + 1  ] = icell ( ix + 1, iy    );
			map_[ imap + 2  ] = icell ( ix + 1, iy + 1);
			map_[ imap + 3  ] = icell ( ix, iy + 1);
			map_[ imap + 4  ] = icell ( ix - 1, iy + 1);
		}
	}
}

void links (int it)
{
	//ZERO head OF CHAiN ARRAY
	for (int icell = 1; icell <= ncell; icell++ )
	{
		head[icell] = -1;
	}
	double celli = M;
	double cell  = BOXL / celli;
	if ( cell < rcut )
		cout << "cell SizE TOO SMALL FOR CUTOFF" << endl;

//SORT ALL ATOMS
	for (int i = 0; i < N; i++)
	{
		int icel =  1 + int ( ( x[i] / BOXL + 0.5 ) * celli )
		            + int ( ( y[i] / BOXL + 0.5 ) * celli ) * M;
		list_[i] = head[icel];
		head[icel] = i;
		if (icel > M * M) cout << "error in icel" << it;
	}
}

void force()
{
	double sigsq  = sigma * sigma;
	double rcutsq = rcut * rcut ;
	double xi, yi, Fxi, Fyi;
	double xij, yij, rijsq, SR2, SR6, Wij, Fij, Fxij, Fyij;
	double r2, r6;

//    ** ZERO FORCES AND POTENTiAL **
	for (int i = 0; i < N ; i++)
	{
		Fx[i] = 0.0;
		Fy[i] = 0.0;
	}
//    ** LOOP OVER ALL cellS **
	for (int icell = 1; icell <= ncell; icell++)
	{
		int  i = head[icell];
//C       ** LOOP OVER ALL MOLECULES iN THE cell **
		while ( i > -1 )
		{
			xi = x[i];
			yi = y[i];
			Fxi = Fx[i];
			Fyi = Fy[i];
//C          ** LOOP OVER ALL MOLECULES BELOW i iN THE CURRENT cell **
			int  j = list_[i];
			while ( j > -1 )
			{
				xij  = x[i] - x[j];
				yij  = y[i] - y[j];

				// periodic boundary conditions
				xij  = xij - BOXL * rint(xij / BOXL);
				yij  = yij - BOXL * rint(yij / BOXL);

				rijsq = xij * xij + yij * yij ;
				if ( rijsq < rcutsq )
				{
					r2 = 1.0 / rijsq;
					r6 = r2 * r2 * r2;
					Wij = (r6 * r6 - 0.5 * r6);
					Fij = Wij * r2;

					Fxij = Fij * xij;
					Fyij = Fij * yij;
					Fxi += Fxij;
					Fyi += Fyij;
					Fx[j] -= Fxij;
					Fy[j] -= Fyij;
				}
				j = list_[j];
			}
			int jcell0, jcell;
			jcell0 = 4 * (icell - 1);
			for ( int nabor = 1; nabor <= 4; nabor++)
			{
				jcell = map_[jcell0 + nabor];
//C             ** LOOP OVER ALL MOLECULES iN NEiGHBOURiNG cellS **
				j = head[jcell];
				while ( j > -1 )
				{
					xij  = x[i] - x[j];
					yij  = y[i] - y[j];

					xij  = xij - BOXL * rint( xij / BOXL);
					yij  = yij - BOXL * rint( yij / BOXL);
					rijsq = xij * xij + yij * yij;
					if ( rijsq < rcutsq )
					{
						r2 = 1.0 / rijsq;
						r6 = r2 * r2 * r2;
						Wij = (r6 * r6 - 0.5 * r6);
						Fij = Wij * r2;

						Fxij = Fij * xij;
						Fyij = Fij * yij;
						Fxi += Fxij;
						Fyi += Fyij;
						Fx[j] -= Fxij;
						Fy[j] -= Fyij;
					}
					j = list_[j];
				}
			}
			Fx[i] = Fxi;
			Fy[i] = Fyi;

			i = list_[i];
		}
	}
	for (int i = 0; i < N; i++)
	{
		Fx[i] = 8 * alpha * Fx[i];
		Fy[i] = 8 * alpha * Fy[i];
	}
}
double DrConstant = sqrt(2 * Dr * dt);
double DConstant = sqrt(2 * D / dt);

void update_position_angle()
{
	for (int i = 0; i < N; i++)
	{
		vx[i] = (Fx[i] + Pe * cos(theta[i])) + DConstant * gau_dist(generator);
		vy[i] = (Fy[i] + Pe * sin(theta[i])) + DConstant * gau_dist(generator);
		x[i] += vx[i] * dt;
		y[i] += vy[i] * dt;
		//periodic boundary condition
		x[i] = x[i] - BOXL * rint(x[i] / BOXL);
		y[i] = y[i] - BOXL * rint(y[i] / BOXL);
		theta[i] +=  DrConstant * gau_dist(generator);
	}
}

void write_VMD_data(ostream& os)
{
	os << N << endl;
	os << "LM" << endl;
	os.setf(ios::fixed, ios::floatfield);
	os.precision(5);
	for (int i = 0; i < N; i++)
	{
		os << "s " << x[i] << " " << y[i] << " " << theta[i] << endl;
	}
	os.flush();
}

int main(int argc, char *argv[])
{
	ofstream o("vmd_data.xyz");
	ofstream out("data.xyz");

	out << "  no of particles::" << N << "  separation: a:: " << a << "  rcut: " << rcut << flush << endl;
	out << "boxlength::" << BOXL << "  M: " << M << "  cell size: L / M::" << L / M << "  cell size: BOXL / M::" << BOXL / M << flush << endl;
	out << " sigma: " << sigma << " MAXIT: " << MAXIT << flush << endl;
	out << "Pe: " << Pe << " phi: " << phi << " alpha: " << alpha << flush << endl;

	double frame = par_maxFrame;
	int tn = MAXIT / frame;

	initialize_square_lattice();
	maps();
	for (int it = 0; it < MAXIT; it++)
	{
		links(it);
		force();
		if (it % tn == 0)
			write_VMD_data(o);
		update_position_angle();
	}
	out << "Time elapsed: " << ( (double(clock()) - start) / CLOCKS_PER_SEC) << endl;
}
