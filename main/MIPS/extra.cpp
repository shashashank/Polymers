#include "extra.h"

// Reads contents of file "file" to array fileLines and returns its length
int readFile(std::string file, std::vector<std::string>& data){
    std::ifstream infile(file);
    std::string line;
    while (getline(infile, line)){
        data.push_back(line);
    }
    return data.size();
}

// Imports to supplied vectors X and Y the final (R=0)
//  or random configuration from data file "oldFile"
void extractConfig(std::string oldFile, std::vector<double>& X, std::vector<double>& Y, std::vector<double>& THETA){
    int N, tPos, fileLength, iters;
    std::vector<std::string> fileLines;
    
    fileLength = readFile(oldFile, fileLines);

    std::istringstream iss(fileLines[0]);
    iss >> N;
    printf("Total number of particles: %i \n", N);
    iters = fileLength / (N + 2);
    tPos = (iters - 1) * (N + 2);
    char s;
    double a, b, c;
    int count = 0;
    if (X.size()==N){
        for (int i = tPos + 2; i < (tPos + 2 + N); ++i){
            std::istringstream iss(fileLines[i]);
            iss >> s >> a >> b;
            X[count] = a;
            Y[count] = b;
            THETA[count] = c;
            ++count;
        }
    }else{
        for (int i = tPos + 2; i < (tPos + 2 + N); ++i){
            std::istringstream iss(fileLines[i]);
            iss >> s >> a >> b;
            X.push_back(a);
            Y.push_back(b);
            THETA.push_back(c);
            ++count;
        }
    }
    printf("Number of particles initialised: %i \n", count);
}