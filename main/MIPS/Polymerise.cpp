/* This compiles to give an executable which you need to run with the following
runtime arguments:
    1) the vmd data file
    2) the number of the particle from the vmd file
    3) the angle
    4) the length of the polymer to be created */

#include "extra.h"		// headers go in here


int N, iters;
std::vector<double> X, Y, THETA;
std::vector<int> polymer;
std::list<int> ascPoly;

double distance(double x1, double y1, double x2, double y2);
int buildPolymer(int start, double angle);
void writePolyVMD(std::ostream& poly);

int main(int argc, char *argv[]){
    std::string vmdData;
    int N, pos;
    double rad;
    std::ofstream polyEle("polymerElements.xyz");

    if (argc<5){
        std::cout << "Please provide the initial monomer particle, and orientation." << std::endl;
        return(1);
    }else{
        std::stringstream(argv[1]) >> vmdData;
        std::stringstream(argv[2]) >> pos;
        std::stringstream(argv[3]) >> rad;
        std::stringstream(argv[4]) >> N;
        rad = rad*M_PI/180;
        printf("Starting Particle:  %i \nAngle: %f \nNo. of Monomers:   %i\n", pos, rad, N);
    }

    extractConfig(vmdData, X, Y, THETA);
    polymer.push_back(pos);
    
    for (size_t i = 0; i < N; i++)
    {
        pos = buildPolymer(pos, rad);
        polymer.push_back(pos);
    }

    writePolyVMD(polyEle);
}

int buildPolymer(int start, double angle)
{
    int nextMonomer;
    double x, y, x1, y1, partDist, partAngle, shortestDist;
    x1 = X[start];
    y1 = Y[start];
    shortestDist = 2;
    for (int i = 0; i < X.size(); ++i){
        x = X[i] - x1;
	    y = Y[i] - y1;
	    partDist = sqrt(x*x + y*y);
        partAngle = atan2(y, x);
        if (partDist < shortestDist && (angle+0.6 > partAngle && angle-0.6 < partAngle)){
            shortestDist = partDist;
            //printf("Start: %i \nDist:  %f \nAngle: %f \ni:    %i\n", start, partDist, partAngle, i);
            nextMonomer = i;
        }
    }
    return nextMonomer;
}

double distance(double x1, double y1, double x2, double y2)
{
	double x = x2 - x1;
	double y = y2 - y1;
	return sqrt(x*x + y*y), atan2(y, x);
}

void writePolyVMD(std::ostream& poly){
    int j = 0;
	poly << X.size() << std::endl;
	poly << "LM" << std::endl;
	poly.setf(std::ios::fixed, std::ios::floatfield);
	poly.precision(5);

    if (ascPoly.size()==0)
    {
        for (int i = 0; i < polymer.size(); i++)
        {
            ascPoly.push_back(polymer[i]);
        }
        ascPoly.sort();
    }

    std::list<int>::iterator ptr = ascPoly.begin();
    for (int i = 0; i < X.size(); i++)
    {
        if (*ptr==i)
        {
            poly << "p "<< X[i] << " " << Y[i] << " " << THETA[i] << std::endl;
            std::advance(ptr, 1);
        }else{
            poly << "s " << X[i] << " " << Y[i] << " " << THETA[i] << std::endl;
        }
    }
    printf("Succesfully written to file.\n");
    
}