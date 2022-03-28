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
void extractConfig(std::string oldFile, std::vector<double>& X, std::vector<double>& Y, int R){
    int N, t_pos, fileLength, iters;
    std::vector<std::string> fileLines;
    
    fileLength = readFile(oldFile, fileLines);

    std::istringstream iss(fileLines[0]);
    iss >> N;
    printf("Total number of particles: %i \n", N);
    iters = fileLength / (N + 2);
    if (R==0){
        t_pos = (iters - 1) * (N + 2);
    }else{
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> distr(1, iters);
        t_pos = distr(gen) * (N + 2);
    }
    char s;
    double a, b;
    int count = 0;
    for (int i = t_pos + 2; i < (t_pos + 2 + N); ++i){
        std::istringstream iss(fileLines[i]);
        iss >> s >> a >> b;
        X[count] = a;
        Y[count] = b;
        ++count;
    }
    printf("Number of particles initialised: %i \n", count);
}

