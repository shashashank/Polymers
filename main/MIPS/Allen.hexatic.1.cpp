#include "hexaticMain.h"
#include "extra.h"
#include "polymer.h"
#include "param.c"

int N1 = 125220;//par_N;
double tt = 0.0;
double phi = par_phi;
double Pe = par_Pe;
double alpha = 6;
constexpr int sigma = 1;
constexpr double Gamma = 1.0;
constexpr double c = 1.0;

constexpr double D = 1;
double Dr = 3 * D / (sigma*sigma);

constexpr double pi = M_PI;
constexpr double TAU = 1;
constexpr double dt = 0.00001 * TAU;
double MAXIT =  par_maxIT * TAU / dt;

/******calculate parameter for cell list*******/
double rcut=sigma*pow(2,1.0/alpha);
int n;//=sqrt(N1);
int N=125220;//n * n;
double L = sqrt(N*pi*sigma*sigma*0.25/phi);
int M = int(L/rcut);
double a=2*L/(2*n+1);
double DrConstant = sqrt(2 * Dr * dt);
double DConstant = sqrt(2 * D / dt);
double BOXL=M*rcut;
int ncell = M*M;
int mapsiz = 4*ncell;
double A = BOXL*BOXL;
/***********  ***********/

std::vector<int> list_(N);
std::vector<double> Fx(N), Fy(N);
std::vector<int> head(ncell + 1);
std::vector<int> map_(mapsiz + 2);
int ix, iy, iz, imap;

std::vector<double> theta(N);
std::vector<double> x(N), y(N);
std::vector<double> vx(N);
std::vector<double> vy(N);
std::vector<int> polyList;
std::list<int> ascendingPolyList;

std::mt19937 generator(time(NULL));
std::normal_distribution<double> gau_dist(0.0, 1.0);
std::uniform_real_distribution<double> uni_dist(-pi / 4, pi / 4);


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

			theta[nptl] = 0;//sqrt(2 * Dr * dt) * gau_dist(generator);
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
			map_[ imap + 3  ] = icell ( ix	  , iy + 1);
			map_[ imap + 4  ] = icell ( ix - 1, iy + 1);
		}
	}
}

void links (int it)
{
	//ZERO head OF CHAiN ARRAY
	for(int icell=1; icell <= ncell; icell++ )
	{
		head[icell] = -1;
	}
	double celli = M;
	double cell  = BOXL / celli;
	if( cell<rcut )
		std::cout<< "cell SizE TOO SMALL FOR CUTOFF\n"<<std::endl;

//SORT ALL ATOMS
	for(int i=0; i<N; i++)
	{
		int icel =  1+int ( ( x[i]/BOXL + 0.5 ) * celli )
		            + int ( ( y[i]/BOXL + 0.5 ) * celli ) * M;
		list_[i] = head[icel];
		head[icel] = i;
		if(icel>M*M) std::cout<<"error in icel\n"<<it<<std::endl;
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
	for(int i =0; iN ; i++)
	{
		Fx[i] = 0.0;
		Fy[i] = 0.0;
	}
//    ** LOOP OVER ALL cellS **
	for(int icell=1; icell<=ncell; icell++)
	{
		int  i = head[icell];
//C       ** LOOP OVER ALL MOLECULES iN THE cell **
		while( i >-1 )
		{
			xi = x[i];
			yi = y[i];
			Fxi = Fx[i];
			Fyi = Fy[i];
//C          ** LOOP OVER ALL MOLECULES BELOW i iN THE CURRENT cell **
			int  j = list_[i];
			while( j >-1 )
			{
				xij  = x[i] - x[j];
				yij  = y[i] - y[j];

				// periodic boundary conditions
				xij  = xij - L * rint(xij / L);
				yij  = yij - L * rint(yij / L);

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

	for (size_t i = 0; i < polyList.size(); i++)
	{
		if (i!=0){
			springForce(x, y, Fx, Fy, polyList[i], polyList[i-1], 1.0, 100); //spring force from previous particle in chain
	  	}
	  	if (i!=(N-1)){
			springForce(x, y, Fx, Fy, polyList[i], polyList[i+1], 1.0, 100); //spring force from next particle in chain
		}
		
	}
	
}

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

void write_VMD_data(std::ostream& os)
{
	os << N << std::endl;
	os << "LM" << std::endl;
	os.setf(std::ios::fixed, std::ios::floatfield);
	os.precision(5);
	for (int i = 0; i < N; i++)
	{
		os << "s " << x[i] << " " << y[i] << " " << theta[i] << std::endl;
	}
	os.flush();
}

int main(int argc, char *argv[])
{
	clock_t start = clock();
	std::ofstream o("vmd_data.xyz");
	std::ofstream out("data.xyz");

	out << "  no of particles::" << N << "  separation: a:: " << a << "  rcut: " << rcut << std::endl;
	out << "boxlength::" << BOXL << "  M: " << M << "  cell size: L / M::" << L / M << "  cell size: BOXL / M::" << BOXL / M << std::endl;
	out << " sigma: " << sigma << " MAXIT: " << MAXIT << std::endl;
	out << "Pe: " << Pe << " phi: " << phi << " alpha: " << alpha << std::endl;

	double frame = par_maxFrame;
	int tn = MAXIT / frame;
	
    printf("Number of particles: %i \n", N);
	if (std::experimental::filesystem::v1::exists("vmd_data_old.xyz")){
		extractConfig("vmd_data_old.xyz", x, y, theta);
		// makePolymer("vmd_data_old.xyz", 51389, 0.45, 50, o, x, y, theta, polyList, ascendingPolyList);
	}
	else initialize_square_lattice();
	maps();

	for (int it = 0; it < MAXIT; it++)
	{
		links(it);
		force();

		if (it % tn == 0){
			write_VMD_data(o);
			//writePolyVMD(o, x, y, theta, polyList, ascendingPolyList);
		}
		update_position_angle();
	}
	out << "Time elapsed: " << ( (double(clock()) - start) / CLOCKS_PER_SEC) << std::endl;
}
