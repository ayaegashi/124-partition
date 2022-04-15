#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <string>
#include <cstdlib>
#include <fstream>
#include <iostream>

using namespace std;

int main(int argc, char *argv[]) {

    //checks that there are command line args
    if (argc != 4)
    {
        cout << "Usage: ./partition flag algorithm inputfile\n";
        return 1;
    }

    // ALGORITHM CODE
    // 0: Karmarkar-Karp
    // 1: Repeated Random
    // 2: Hill Climbing
    // 3: Simulated Annealing
    // 11: Prepartitioned Repeated Random
    // 12: Prepartitioned Hill Climbing
    // 13: Prepartitioned Simulated Annealing

    int flag = (int) strtol(argv[1], NULL, 10);
    int alg = (int) strtol(argv[2], NULL, 10);
    char *inputfile = argv[3];



    return 0;
};

// Karmarkar-Karp Algorithm
int* kk(int* seq){
    int* signs;

    return signs
};