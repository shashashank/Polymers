#include "extra.h"		// headers go in here
#include "poly_param.c" // parameters for simulation

double dt = 1e-3;	            // time evolution step size
double l = 1.0;                 // length of monomer
double intRad = l;            // interaction force cutoff radius
int k = 100;                    //
double sqrt_dt = sqrt(dt);
double MAXIT =  parMaxIT / dt;  // maximum number of iterations for the simulation
double start_time;              // stores time of start of execution
double frame = parMaxFrame;     // number of times to save data
int tn = MAXIT / frame;

std::mt19937 generator(time(NULL));
std::normal_distribution<double> gau_dist(0.0, 1.0);

std::vector<double> Px(N), Py(N);
std::vector<double> Fx(N), Fy(N);
//std::vector<double> intrX(N), intrY(N);
std::vector<double> e2eArray(frame);
// std::vector<omp_lock_t> iXLock(N), iYLock(N);

double spring(double x, double y, double length, int spring_constant);
void animate(void);
void initialize(void);
void write_VMD_data(std::ostream& os, int it);
void noise(void);
double e2e_distance(double x_ini, double x_fin, double y_ini, double y_fin);
void SForce(int P1, int P2);
void IForce(int P1, int P2, double *intX, double *intY);
void writeE2eData(std::ostream& os, int it);

class progress_bar
{
    static const auto overhead = sizeof " [100%]";

    std::ostream& os;
    const std::size_t bar_width;
    std::string message;
    const std::string full_bar;

 public:
    progress_bar(std::ostream& os, std::size_t line_width,
                 std::string message_, const char symbol = '.')
        : os{os},
          bar_width{line_width - overhead},
          message{std::move(message_)},
          full_bar{std::string(bar_width, symbol) + std::string(bar_width, ' ')}
    {
        if (message.size()+1 >= bar_width || message.find('\n') != message.npos) {
            os << message << '\n';
            message.clear();
        } else {
            message += ' ';
        }
        write(0.0);
    }

    // not copyable
    progress_bar(const progress_bar&) = delete;
    progress_bar& operator=(const progress_bar&) = delete;

    ~progress_bar()
    {
        write(1.0);
        os << '\n';
    }

    void write(double fraction);
};

void progress_bar::write(double fraction)
{
    // clamp fraction to valid range [0,1]
    if (fraction < 0)
        fraction = 0;
    else if (fraction > 1)
        fraction = 1;

    auto width = bar_width - message.size();
    auto offset = bar_width - static_cast<unsigned>(width * fraction);

    os << '\r' << message;
    os.write(full_bar.data() + offset, width);
    os << " [" << std::setw(3) << static_cast<int>(100*fraction) << "%] " << std::flush;
}

int main(int argc, char *argv[]){
    start_time = omp_get_wtime();
	std::ofstream out("data_poly.xyz");
	std::ofstream o("vmd_data_poly.xyz");
	std::ofstream e2e("e2e.d");
	std::string old = "old.xyz";

    omp_set_dynamic(1);
	omp_set_num_threads(omp_get_max_threads());

    out << "Number of particles: " << N << std::endl;
    out << "The total number of timesteps: " <<  parMaxIT
        << "/" << dt << " = " << MAXIT << std::endl;

    if (std::experimental::filesystem::v1::exists(old)){
        std::cout << "Importing configuration from older file." << std::endl;
        extractConfig(old, Px, Py, 1);
	}
	else initialize();

	progress_bar progress{std::clog, 70u, "Working"};

// #pragma omp parallel for
// 	for (int ilock = 0; ilock<N; ilock++){
// 		omp_init_lock(&iXLock[ilock]);
// 		omp_init_lock(&iYLock[ilock]);
// 	}
	int f = 0; 

	for (int it = 0; it < MAXIT; it++){
		if ((it % tn == 0) && (writeFlag==1)){
			progress.write(f/frame);
            write_VMD_data(o, it);
			e2e << e2e_distance(Px[0], Px[N-1], Py[0], Py[N-1]) << std::endl;
			f += 1;
		}
		animate();
	}

    // #pragma omp parallel for
	// for (int ilock = 0; ilock<N; ilock++){
	// 	omp_destroy_lock(&iXLock[ilock]);
	// 	omp_destroy_lock(&iYLock[ilock]);
	// }
	
	out << "Time elapsed: " << (omp_get_wtime() - start_time) << std::endl;

}

void initialize(void){
	for (int i=0; i<N; i++)
	{
		Px[i] = double(i * l);
		Py[i] = 0.0;
	}
}

void animate(void){
	double apart;
	double intrXThread[N] = {0};
	double intrYThread[N] = {0};
// add schedule(guided, 10) for longer chains
#pragma omp parallel for schedule(static, 8)
	for (int i=0; i<N; i++)
	{
	  	Fx[i] = 0.0; Fy[i] = 0.0;
	  	if (i!=0){
			SForce(i, i-1); //spring force from previous particle in chain
	  	}
	  	if (i!=(N-1)){
			SForce(i, i+1); //spring force from next particle in chain
		}
    }
#pragma omp parallel for schedule(guided, 8) reduction(+:intrXThread[:N],intrYThread[:N])
    for (int i=0; i<N; i++){
		for (int j=i+2; j<N; j++){
			apart = e2e_distance(Px[i], Px[j], Py[i], Py[j]) + intrFlag;
			if (apart<intRad){
				IForce(i, j, intrXThread, intrYThread);
			}
		}
	}
// add schedule(static, 100) for longer chains
#pragma omp parallel for schedule(static, 8)
	for (int i=0; i<N; i++){
		Px[i] += (Fx[i] + intrXThread[i])*dt + gau_dist(generator)*sqrt_dt;
		Py[i] += (Fy[i] + intrYThread[i])*dt + gau_dist(generator)*sqrt_dt;
		intrXThread[i] = 0.0;
		intrXThread[i] = 0.0;
	}
}
double spring(double x, double y, double length, int spring_constant) {
    double K = double(spring_constant);
    double distance = sqrt(x*x + y*y);
    return K*(distance - length);
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
void IForce(int P1, int P2, double *intX, double *intY){
	double x, y, phi, F, Fx, Fy;
	x = Px[P2] - Px[P1];
	y = Py[P2] - Py[P1];
	F = spring(x, y, intRad, k);
	phi = atan2(y, x);
	Fx = F*cos(phi);
	Fy = F*sin(phi);
	//phi2 = atan2((-1.0*y), (-1.0*x));

	//omp_set_lock(&iXLock[P1]);
	intX[P1] += Fx;
	//omp_unset_lock(&iXLock[P1]);

	//omp_set_lock(&iYLock[P1]);
	intY[P1] += Fy;
	//omp_unset_lock(&iYLock[P1]);

	//omp_set_lock(&iXLock[P2]);
	intX[P2] -= Fx;
	//omp_unset_lock(&iXLock[P2]);

	//omp_set_lock(&iYLock[P2]);
	intY[P2] -= Fy;
	//omp_unset_lock(&iYLock[P2]);
}

void write_VMD_data(std::ostream& os, int it)
{
	os << N << std::endl;
	os << it << std::endl;
	os.setf(std::ios::fixed, std::ios::floatfield);
	os.precision(5);
	for (int i = 0; i < N; i++)
	{
		os << "s " << Px[i] << " " << Py[i] << std::endl;
	}
	os.flush();
}

double e2e_distance(double x_ini, double x_fin, double y_ini, double y_fin)
{
	double X = x_fin - x_ini;
	double Y = y_fin - y_ini;
	return sqrt(X*X + Y*Y);
}
