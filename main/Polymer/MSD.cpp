#include "extra.h"

void meanSqauredDisplacement();
void radiusOfGyration();
void writeData();


int N, iters;
vector<vector<double> > X;
vector<vector<double> > Y;
vector<vector<double> > msdArray;
vector<vector<double> > centreOfMass;
vector<double> gyrationRadii;
ofstream rgData("radiiOfGyration.d");
ofstream msdData("msdVals.d");

int main(int argc, char *argv[]){
    string vmdData = "vmd_data_poly.xyz";
    int startIter = 1000;
    int fileLength, tPos, count, lenCheck;
    vector<string> fileLines;

    fileLength = readFile(vmdData, fileLines);

    istringstream iss(fileLines[0]);
    iss >> N;

    iters = (fileLength / (N + 2)) - startIter;

    // X and Y will hold X and Y positions of all the time steps
    X.resize(iters, vector<double>(N));
    Y.resize(iters, vector<double>(N));
    msdArray.resize(N-1, vector<double>(N));
    centreOfMass.resize(iters, vector<double>(2));
    gyrationRadii.resize(iters);

    char s;
    tPos = 0;
    lenCheck = startIter*(2 + N);
    double a, b;
    while(tPos < iters){
        count = 0;
        for (int i = lenCheck + 2; i < (lenCheck + 2 + N); ++i){
            istringstream iss(fileLines[i]);
            iss >> s >> a >> b;
            X[tPos][count] = a;
            Y[tPos][count] = b;
            ++count;
        }
        ++tPos;
        lenCheck = (tPos + startIter)*(2 + N);
    }

    for (int i = 0; i<iters; i++){
        centreOfMass[i][0] = 0.0;
        centreOfMass[i][1] = 0.0;
        for (int j = 0; j<N; j++){
            centreOfMass[i][0] = centreOfMass[i][0]+X[i][j];
            centreOfMass[i][1] = centreOfMass[i][1]+Y[i][j];
        }
        centreOfMass[i][0] = centreOfMass[i][0]/N;
        centreOfMass[i][1] = centreOfMass[i][1]/N;
    }

    meanSqauredDisplacement();
    radiusOfGyration();
    writeData();

}

// Calculates the MSD(tau), with increasing tau
void meanSqauredDisplacement(){
    for (int i = 0; i<(N-1); i++){
        for (int j = 0; j<N; j++){
            msdArray[i][j] = 0.0;
            for (int k = 0; k<iters; k++){
                if ((i+k)<N){
                    msdArray[i][j] = msdArray[i][j] + (X[k+i][j] - X[k][j])*(X[k+i][j] - X[k][j]);
                    msdArray[i][j] = msdArray[i][j] + (Y[k+i][j] - Y[k][j])*(Y[k+i][j] - Y[k][j]);
                }
            }
            msdArray[i][j] = msdArray[i][j] / (N - i + 1);
        }
    }
}

// Calculates radius of gyration
void radiusOfGyration(){
    for (int i = 0; i < iters; i++){
        gyrationRadii[i] = 0.0;
        for(int j = 0; j<N; j++){
            gyrationRadii[i] = (X[i][j] - centreOfMass[i][0])*(X[i][j] - centreOfMass[i][0]);
            gyrationRadii[i] = (Y[i][j] - centreOfMass[i][1])*(Y[i][j] - centreOfMass[i][1]);
        }
        gyrationRadii[i] = gyrationRadii[i]/N;
    }
}

void writeData(){
	rgData << N << endl;
	rgData.setf(ios::fixed, ios::floatfield);
	rgData.precision(5);
    for (int i=0; i<iters; i++){
        rgData << gyrationRadii[i] << endl;
    }
	rgData.flush();

   	msdData << N << endl;
	msdData.setf(ios::fixed, ios::floatfield);
	msdData.precision(5);
    for (int i=0; i<N-1; i++){
        for (int j = 0; j<N; j++){
            msdData << msdArray[i][j] << "\t";
        }
        msdData << endl;
    }
	msdData.flush();
}