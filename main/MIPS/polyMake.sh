#!/bin/zsh
clear
startN=460070
angle=150
lenN=200
cd /home/shashank/Cross/Polymers/main/MIPS
g++ -std=c++17 -O -fopenmp extra.cpp Polymerise.cpp -lstdc++fs -o polymerCreator
cp polymerCreator ../../../PolymerChainCreator
cd /home/shashank/Cross/PolymerChainCreator
./polymerCreator vmd_data.xyz $startN $angle $lenN
