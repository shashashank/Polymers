#include "extra.h"


// Reads contents of file "file" to array fileLines and returns its length
int readFile(string file, vector<string>& data){
    ifstream infile(file);
    string line;
    while (getline(infile, line)){
        data.push_back(line);
    }
    return data.size();
}

// Imports to supplied vectors X and Y the final (R=0)
//  or random configuration from data file "oldFile"
void extractConfig(string oldFile, vector<double>& X, vector<double>& Y, int R){
    int N, t_pos, fileLength, iters;
    vector<string> fileLines;
    
    fileLength = readFile(oldFile, fileLines);

    istringstream iss(fileLines[0]);
    iss >> N;
    printf("Total number of particles: %i \n", N);
    iters = fileLength / (N + 2);
    if (R==0){
        t_pos = (iters - 1) * (N + 2);
    }else{
        random_device rd;
        mt19937 gen(rd());
        uniform_int_distribution<> distr(1, iters);
        t_pos = distr(gen) * (N + 2);
    }
    char s;
    double a, b;
    int count = 0;
    for (int i = t_pos + 2; i < (t_pos + 2 + N); ++i){
        istringstream iss(fileLines[i]);
        iss >> s >> a >> b;
        X[count] = a;
        Y[count] = b;
        ++count;
    }
    printf("Number of particles initialised: %i \n", count);
}

