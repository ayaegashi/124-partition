#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <string>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
#include <tuple>
#include <algorithm>
#include <random>

using namespace std;

random_device rd;
mt19937 gen(rd());

const int MAX_ITER = 100;

// Algorithms
int kk(vector<uint64_t> seq);
int repeatedRandom(vector<uint64_t> A, int n, bool isSequence);
int hillClimbing(vector<uint64_t> A, int n, bool isSequence);
int simulatedAnnealing(vector<uint64_t> A, int n, bool isSequence);
int T(int i);

// For Prepartitioning
int calcResiduePrepartitioning(vector<uint64_t> A, vector<int> S);
vector<int> generateRandomPrepartitioningSoln(int n);
vector<int> generateRandomPrepartitioningMove(vector<int> prev, int n);

// Generating Instances
vector<uint64_t> generateRandomInstance(int n);
int generateRandomInt(int low, int high);

// THINGS JENNS NOT SURE ABOUT: should the residue be an unit64_t as well??? instead of just an int???
// I'm also not sure if I implemented 64 bit ints correctly lol

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

    vector<uint64_t> sequence;
    // Construct vector from inputfile
    if (flag == 0) {
        ifstream file (inputfile, ios::in);
        if (file.fail()) {
            printf("File could not be opened.\n");
            return 0;
        }
        uint64_t n;
        while (file >> n){
            sequence.push_back(n);
        }

        file.close();
    }
    // Generate x random instances from function
    else if (flag == 1) {
        for (int i = 0; i < 1; i++) {
            sequence = generateRandomInstance(100);
            int n = sequence.size();

            int kkResidue = kk(sequence);
            printf("kk residue: %i\n", kkResidue);

            int repeatedRandomResidue = repeatedRandom(sequence, n, false);
            printf("repeated random residue: %i\n", repeatedRandomResidue);

            int hillClimbingResidue = hillClimbing(sequence, n, false);
            printf("hill climbing residue: %i\n", hillClimbingResidue);

            int simulatedAnealingResidue = simulatedAnnealing(sequence, n, false);
            printf("simulated annealing residue: %i\n", simulatedAnealingResidue);
        }
    }

    int n = sequence.size();

    // for (int i = 0; i < n; i ++) {
    //     printf("%llu ", sequence[i]);
    // }

    if (flag == 0 && alg == 0){
        int difference = kk(sequence);
        printf("\n%i\n", difference);
    }
    else if (flag == 0 && alg == 11) {
        int residue = repeatedRandom(sequence, n, false);
        printf("\n %i\n", residue);
    }
    else if (flag == 0 && alg == 12) {
        int residue = hillClimbing(sequence, n, false);
        printf("\n %i\n", residue);
    }
    else if (flag == 0 && alg == 13) {
        int residue = simulatedAnnealing(sequence, n, false);
        printf("\n %i\n", residue);
    }

    return 0;
};

// IGNORE FOR NOW DOESNT WORK!
// Karmarkar-Karp Algorithm --> returns just difference
int kk(vector<uint64_t> seq){
    make_heap(seq.begin(), seq.end());

    int largest = 0;
    int nextLargest = 0;
    do {
        int largest = seq.front();
        pop_heap (seq.begin(),seq.end());  // puts the largest element at the end of the vector
        seq.pop_back();  // removes the last element of the vector (what we just popped from the heap)
        int nextLargest = seq.front();
        pop_heap (seq.begin(),seq.end());  // puts the largest element at the end of the vector
        seq.pop_back();  // removes the last element of the vector (what we just popped from the heap)

        int difference = largest - nextLargest;
        if (difference > 0){
            seq.push_back(difference);  // add difference to end of the vector
            push_heap (seq.begin(), seq.end());
        }

    } while (seq.size() > 1);

    if (seq.size() == 1){
        return seq.front();
    } else {
        return 0;
    }

    // This doesn't actually do anything. Just was getting an annoying compiler warning without this.
    largest += 1;
    nextLargest += 1;
};

// Returns best residue from repeated random algorithm
int repeatedRandom(vector<uint64_t> A, int n, bool isSequence) {
    if (isSequence) {
        //TODO: For Ayana: sequence code here
        return -1;
    // For Prepartitioning representation
    } else {
        vector<int> bestSoln = generateRandomPrepartitioningSoln(n);
        int bestResidue = calcResiduePrepartitioning(A, bestSoln);
        for (int i = 0; i < MAX_ITER; i++) {
            vector<int> soln = generateRandomPrepartitioningSoln(n);
            int residue = calcResiduePrepartitioning(A, soln);
            if (residue < bestResidue) {
                bestSoln = soln;
                bestResidue = residue;
            }
        }
        return bestResidue;
    }
}

// Returns best residue from hill climbing algorithm
int hillClimbing(vector<uint64_t> A, int n, bool isSequence) {
    if (isSequence) {
        //TODO: For Ayana: Sequence Code Here
        return -1;
    } else {
        vector<int> startSoln = generateRandomPrepartitioningSoln(n);
        int bestResidue = calcResiduePrepartitioning(A, startSoln);
        for (int i = 0; i < MAX_ITER; i++) {
            vector<int> neighbor = generateRandomPrepartitioningMove(startSoln, n);
            int residue = calcResiduePrepartitioning(A, neighbor);
            if (residue < bestResidue) {
                bestResidue = residue;
                startSoln = neighbor;
            }
        }
        return bestResidue;
    }
}

int simulatedAnnealing(vector<uint64_t> A, int n, bool isSequence) {
    if (isSequence) {
        //TODO: For Ayana: Sequence Code Here
        return -1;
    } else {
        vector<int> currSoln = generateRandomPrepartitioningSoln(n);
        vector<int> bestSoln = currSoln;
        int currResidue = calcResiduePrepartitioning(A, currSoln);
        int bestResidue = currResidue;
        for (int i = 0; i < MAX_ITER; i++) {
            vector<int> neighbor = generateRandomPrepartitioningMove(currSoln, n);
            int residue = calcResiduePrepartitioning(A, neighbor);
            // if a better residue neighbor OR with some probability switch to it anyway
            if (residue < currResidue || (double) rand() / RAND_MAX < exp((uint64_t) (currResidue - residue) / T(i))) {
                currSoln = neighbor;
                currResidue = residue;
            } 

            if (currResidue < bestResidue) {
                bestResidue = currResidue;
                bestSoln = currSoln;
            }

        }
        return bestResidue;
    }
}

int T(int i) {
    int iter = floor(i/300);
    return pow(10, 10) * pow(0.8, iter);
}

// Generates a random preparitioning solution
vector<int> generateRandomPrepartitioningSoln(int n) {
    vector<int> soln;
    for (int i = 0; i < n; i++) {
        soln.push_back(generateRandomInt(1, n));
    }
    return soln;
}

// Given problem instance A and prepartioning solution S, calculates the residue with KK
int calcResiduePrepartitioning(vector<uint64_t> A, vector<int> S) {
    int n = S.size();
    // Initialize A_prime
    vector<uint64_t> A_prime;
    for (int i = 0; i < n; i++) {
        A_prime.push_back(0);
    }

    // Generate A_i from prepartioning S
    for (int i = 0; i < n; i++) {
        int new_idx = S[i];
        A_prime[new_idx] += A[i];
    }

    // Calculate Residue with KK
    return kk(A_prime);
}

// Takes previous solution and randomly generates next solution to explore out of neighbor state space
vector<int> generateRandomPrepartitioningMove(vector<int> prev, int n) {
    int i = generateRandomInt(1, n);
    int j = generateRandomInt(1, n);
    while (prev[i] == j) {
        i = generateRandomInt(1, n);
        j = generateRandomInt(1, n);
    }
    prev[i] = j;
    return prev;
}

// Generates random instance A of partitioning problem
vector<uint64_t> generateRandomInstance(int n) {
    vector<uint64_t> instance;
    for (int i = 0; i < n; i++) {
        instance.push_back((uint64_t) generateRandomInt(1, pow(10, 12)));
    }
    return instance;
}

// Generates random 64 bit integer between low and high inclusive
int generateRandomInt(int low, int high) {
    std::uniform_int_distribution<> dist(low, high);
    return dist(gen);
}