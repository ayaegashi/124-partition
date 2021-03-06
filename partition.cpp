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
#include <chrono>

using namespace std;

// Max Heap data structure
struct MaxHeap {
    vector<uint64_t> A;

    int SIZE() {
        return A.size();
    }

    int PARENT_INDEX(int i) {
        if (i == 0) {
            // root node
            return 0;
        } else {
            return (i - 1) / 2;
        }
    }

    int LEFT_CHILD(int i) {
        return (2*i + 1);
    };

    int RIGHT_CHILD(int i) {
        return (2*i + 2);
    };

    void SWAP(uint64_t& i, uint64_t& j) {
        uint64_t temp = i;
        i = j;
        j = temp;
    }

    void BUBBLE_UP(int i) {
        int p = PARENT_INDEX(i);

        // violation of max-heap property
        if (A[p] < A[i]) {
            SWAP(A[i], A[p]);
            BUBBLE_UP(p);
        }
    }

    void HEAPIFY(int i) {
        int l = LEFT_CHILD(i);
        int r = RIGHT_CHILD(i);

        int largest = i;

        if (l < SIZE() && A[l] > A[i]) {
            largest = l;
        }
        if (r < SIZE() && A[r] > A[largest]) {
            largest = r;
        }

        if (largest != i) {
            SWAP(A[i], A[largest]);
            HEAPIFY(largest);
        }
    }

    uint64_t DELETE_MAX() {
        if (A.size() == 0) {
            printf("Error: heap is empty, cannot delete_max");
            return -1;
        }

        uint64_t maxElt = A[0];

        A[0] = A.back();
        A.pop_back();
        HEAPIFY(0);

        return maxElt;
    };

    void INSERT(uint64_t elt) {
        A.push_back(elt);

        int newEltIndex = SIZE() - 1;
        BUBBLE_UP(newEltIndex);
    };

    void PRINT() {
        printf("\n");
        for (int i = 0; i < SIZE(); i++) {
            printf("%llu ", A[i]);
        }
        printf("\n");
    }
};



random_device rd;
mt19937 gen(rd());

const int MAX_ITER = 25000;

// Algorithms
uint64_t kk(vector<uint64_t> seq);
uint64_t repeatedRandom(vector<uint64_t> A, vector<int> startS, int n, bool isSequence);
uint64_t hillClimbing(vector<uint64_t> A, vector<int> startS, int n, bool isSequence);
uint64_t simulatedAnnealing(vector<uint64_t> A, vector<int> startS, int n, bool isSequence);
double T(int i);

// For Prepartitioning
uint64_t calcResiduePrepartitioning(vector<uint64_t> A, vector<int> S);
vector<int> generateRandomPrepartitioningSoln(int n);
vector<int> generateRandomPrepartitioningMove(vector<int> prev, int n);

// For Sequence Representation
uint64_t calcResidueSequence(vector<uint64_t> A, vector<int> seq);
vector<int> generateRandomSequenceSoln(int n);
vector<int> generateRandomSequenceMove(vector<int> prevSeq);

// Generating Instances
vector<uint64_t> generateRandomInstance(int n);
uint64_t generateRandomInt(uint64_t low, uint64_t high);

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
        uint64_t num;
        while (file >> num){
            sequence.push_back(num);
        }

        file.close();

        int n = sequence.size();

        if (alg == 0){
            uint64_t difference = kk(sequence);
            printf("\n%llu\n", difference);
        }
        else if (alg == 1) {
            vector<int> startS = generateRandomSequenceSoln(n);
            uint64_t residue = repeatedRandom(sequence, startS, n, true);
            printf("\n %llu\n", residue);
        }
        else if (alg == 2) {
            vector<int> startS = generateRandomSequenceSoln(n);
            uint64_t residue = hillClimbing(sequence, startS, n, true);
            printf("\n %llu\n", residue);
        }
        else if (alg == 3) {
            vector<int> startS = generateRandomSequenceSoln(n);
            uint64_t residue = simulatedAnnealing(sequence, startS, n, true);
            printf("\n %llu\n", residue);
        }
        else if (alg == 11) {
            vector<int> startS = generateRandomPrepartitioningSoln(n);
            uint64_t residue = repeatedRandom(sequence, startS, n, false);
            printf("\n %llu\n", residue);
        }
        else if (alg == 12) {
            vector<int> startS = generateRandomPrepartitioningSoln(n);
            uint64_t residue = hillClimbing(sequence, startS, n, false);
            printf("\n %llu\n", residue);
        }
        else if (alg == 13) {
            vector<int> startS = generateRandomPrepartitioningSoln(n);
            uint64_t residue = simulatedAnnealing(sequence, startS, n, false);
            printf("\n %llu\n", residue);
        }

    }
    // Generate x random instances from function
    else if (flag == 1) {
        using std::chrono::high_resolution_clock;
        using std::chrono::duration;
        using std::chrono::milliseconds;
        // Open output file
        ofstream myfile;
        myfile.open("data 3.txt", ofstream::app);
        myfile << "trial, Karmarkar-Karp, SEQ Start Random Residue, SEQ Repeated Random, SEQ Hill Climbing, SEQ Sim. Annealing, PP Start Random, PP Repeated Random, PP Hill Climbing, PP Sim. Annealing, KK Time, SRR Time, SHC Time, SSA Time, PRR Time, PHC Time, PSA Time" << endl;

        for (int i = 0; i < 50; i++) {
            myfile << i << ", ";
            sequence = generateRandomInstance(100);
            int n = sequence.size();

            std::chrono::high_resolution_clock::time_point start = high_resolution_clock::now();
            uint64_t kkResidue = kk(sequence);
            std::chrono::high_resolution_clock::time_point end = high_resolution_clock::now();

            duration<double, milli> kkTime = end - start;
            printf("kk residue: %llu\n", kkResidue);
            myfile << kkResidue << ", ";

            vector<int> startSequenceSol = generateRandomSequenceSoln(n);
            vector<int> startPrepartSol = generateRandomPrepartitioningSoln(n);

            // Sequencing 

            // Start Random
            uint64_t startRandomSeqResidue = calcResidueSequence(sequence, startSequenceSol);
            printf("SEQUENCE start random residue: %llu\n", startRandomSeqResidue);
            myfile << startRandomSeqResidue << ", ";

            // RR
            start = high_resolution_clock::now();
            uint64_t repeatedRandomSeqResidue = repeatedRandom(sequence, startSequenceSol, n, true);
            end = high_resolution_clock::now();

            duration<double, milli> SRRTime = end - start;
            printf("SEQUENCE repeated random residue: %llu\n", repeatedRandomSeqResidue);
            myfile << repeatedRandomSeqResidue << ", ";

            // HC
            start = high_resolution_clock::now();
            uint64_t hillClimbingSeqResidue = hillClimbing(sequence, startSequenceSol, n, true);
            end = high_resolution_clock::now();

            duration<double, milli> SHCTime = end - start;
            printf("SEQUENCE hill climbing residue: %llu\n", hillClimbingSeqResidue);
            myfile << hillClimbingSeqResidue << ", ";

            // SA
            start = high_resolution_clock::now();
            uint64_t simulatedAnealingSeqResidue = simulatedAnnealing(sequence, startSequenceSol, n, true);
            end = high_resolution_clock::now();

            duration<double, milli> SSATime = end - start;
            printf("SEQUENCE simulated annealing residue: %llu\n", simulatedAnealingSeqResidue);
            myfile << simulatedAnealingSeqResidue << ", ";

            // Prepartitioning

            // Start Random
            uint64_t startPrepartResidue = calcResiduePrepartitioning(sequence, startPrepartSol);
            printf("PREPARTITION start random residue: %llu\n", startPrepartResidue);
            myfile << startPrepartResidue << ", ";

            // RR
            start = high_resolution_clock::now();
            uint64_t repeatedRandomPrepartResidue = repeatedRandom(sequence, startPrepartSol, n, false);
            end = high_resolution_clock::now();

            duration<double, milli> PRRTime = end - start;
            printf("PREPARTITION repeated random residue: %llu\n", repeatedRandomPrepartResidue);
            myfile << repeatedRandomPrepartResidue << ", ";

            // HC
            start = high_resolution_clock::now();
            uint64_t hillClimbingPrepartResidue = hillClimbing(sequence, startPrepartSol, n, false);
            end = high_resolution_clock::now();

            duration<double, milli> PHCTime = end - start;
            printf("PREPARTITION hill climbing residue: %llu\n", hillClimbingPrepartResidue);
            myfile << hillClimbingPrepartResidue << ", ";

            // SA
            start = high_resolution_clock::now();
            uint64_t simulatedAnealingPrepartResidue = simulatedAnnealing(sequence, startPrepartSol, n, false);
            end = high_resolution_clock::now();

            duration<double, milli> PSATime = end - start;
            printf("PREPARTITION simulated annealing residue: %llu\n", simulatedAnealingPrepartResidue);
            myfile << simulatedAnealingPrepartResidue << ", ";

            // Printing Runtimes
            myfile << kkTime.count() << ",  ";
            myfile << SRRTime.count() << ",  ";
            myfile << SHCTime.count() << ",  ";
            myfile << SSATime.count() << ",  ";
            myfile << PRRTime.count() << ",  ";
            myfile << PHCTime.count() << ",  ";
            myfile << PSATime.count() << endl;
        
        }
        myfile.close();
    }

    return 0;
};

// Karmarkar-Karp Algorithm --> returns just difference
uint64_t kk(vector<uint64_t> seq){
    // Make heap
    int n = seq.size();
    MaxHeap seqHeap;
    
    for (int i = 0; i < n; i++) {
        seqHeap.INSERT(seq[i]);
    }

    uint64_t largest = (uint64_t) 0;
    uint64_t nextLargest = (uint64_t) 0;
    do {
        uint64_t largest = seqHeap.DELETE_MAX();
        uint64_t nextLargest = seqHeap.DELETE_MAX();

        uint64_t difference = largest - nextLargest;
        if (difference > 0){
            seqHeap.INSERT(difference);  // add difference to the heap
        }

    } while (seqHeap.SIZE() > 1);

    if (seqHeap.SIZE() == 1){
        return seqHeap.DELETE_MAX();
    } else {
        return 0;
    }

    // This doesn't actually do anything. Just was getting an annoying compiler warning without this.
    largest += 1;
    nextLargest += 1;
};

// Returns best residue from repeated random algorithm
uint64_t repeatedRandom(vector<uint64_t> A, vector<int> startS, int n, bool isSequence) {
    // Sequence of +1 and -1 representation
    if (isSequence) {
        vector<int> bestSoln = startS;
        uint64_t bestRes = calcResidueSequence(A, bestSoln);
        for (int i = 0; i < MAX_ITER; i++) {
            vector<int> newSoln = generateRandomSequenceSoln(n); // Completely random (non-neighbor) solution
            uint64_t newResidue = calcResidueSequence(A, newSoln);
            if (newResidue < bestRes) {
                bestRes = newResidue;
                bestSoln = newSoln;
            }
        }
        return bestRes;
    // For Prepartitioning representation
    } else {
        vector<int> bestSoln = startS;
        uint64_t bestResidue = calcResiduePrepartitioning(A, bestSoln);
        for (int i = 0; i < MAX_ITER; i++) {
            vector<int> soln = generateRandomPrepartitioningSoln(n);
            uint64_t residue = calcResiduePrepartitioning(A, soln);
            if (residue < bestResidue) {
                bestSoln = soln;
                bestResidue = residue;
            }
        }
        return bestResidue;
    }
}

// Returns best residue from hill climbing algorithm
uint64_t hillClimbing(vector<uint64_t> A, vector<int> startS, int n, bool isSequence) {
    if (isSequence) {
        vector<int> bestSoln = startS;
        uint64_t bestRes = calcResidueSequence(A, bestSoln);
        for (int i = 0; i < MAX_ITER; i++) {
            vector<int> neighbor = generateRandomSequenceMove(bestSoln);
            uint64_t neighborRes = calcResidueSequence(A, neighbor);
            if (neighborRes < bestRes) {
                bestRes = neighborRes;
                bestSoln = neighbor;
            }
        }
        return bestRes;
    } else {
        vector<int> currSoln = startS;
        uint64_t bestResidue = calcResiduePrepartitioning(A, currSoln);
        for (int i = 0; i < MAX_ITER; i++) {
            vector<int> neighbor = generateRandomPrepartitioningMove(currSoln, n);
            uint64_t residue = calcResiduePrepartitioning(A, neighbor);
            if (residue < bestResidue) {
                bestResidue = residue;
                currSoln = neighbor;
            }
        }
        return bestResidue;
    }
}

// Returns best residue from simulated annealing algorithm
uint64_t simulatedAnnealing(vector<uint64_t> A, vector<int> startS, int n, bool isSequence) {
    if (isSequence) {
        vector<int> currSoln = startS;
        vector<int> bestSoln = startS;
        uint64_t currResidue = calcResidueSequence(A, currSoln);
        uint64_t bestResidue = currResidue;
        for (int i = 0; i < MAX_ITER; i++) {
            vector<int> neighbor = generateRandomSequenceMove(currSoln);
            uint64_t residue = calcResidueSequence(A, neighbor);
            // if a better residue neighbor OR with some probability switch to it anyway
            if (residue < currResidue || (double) rand() / RAND_MAX < exp((uint64_t) (currResidue - residue) / T(i))) {
                currSoln = neighbor;
                currResidue = residue;
            } 
            // check if current is better than the best, then keep track of best
            if (currResidue < bestResidue) {
                bestResidue = currResidue;
                bestSoln = currSoln;
            }
        }
        return bestResidue;
    } else {
        vector<int> currSoln = startS;
        vector<int> bestSoln = startS;
        uint64_t currResidue = calcResiduePrepartitioning(A, currSoln);
        uint64_t bestResidue = currResidue;
        for (int i = 0; i < MAX_ITER; i++) {
            vector<int> neighbor = generateRandomPrepartitioningMove(currSoln, n);
            uint64_t residue = calcResiduePrepartitioning(A, neighbor);
            // if a better residue neighbor OR with some probability switch to it anyway
            if (residue < currResidue || (double) rand() / RAND_MAX < exp((uint64_t) (currResidue - residue) / T(i))) {
                currSoln = neighbor;
                currResidue = residue;
            } 
            // check if current is better than the best, then keep track of best
            if (currResidue < bestResidue) {
                bestResidue = currResidue;
                bestSoln = currSoln;
            }
        }
        return bestResidue;
    }
}

// Cooling function for simulated annealing
double T(int i) {
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
uint64_t calcResiduePrepartitioning(vector<uint64_t> A, vector<int> S) {
    int n = S.size();
    // Initialize A_prime
    vector<uint64_t> A_prime;
    for (int i = 0; i < n; i++) {
        A_prime.push_back((uint64_t) 0);
    }

    // Generate A_i from prepartioning S
    for (int i = 0; i < n; i++) {
        int new_idx = S.at(i);
        // printf("%i\n", new_idx);
        A_prime[new_idx-1] += A.at(i);
    }

    // Calculate Residue with KK
    return kk(A_prime);
}

// Takes previous solution and randomly generates next solution to explore out of neighbor state space
vector<int> generateRandomPrepartitioningMove(vector<int> prev, int n) {
    int i = generateRandomInt(0, n-1);
    int j = generateRandomInt(1, n);
    while (prev.at(i) == j) {
        i = generateRandomInt(0, n-1);
        j = generateRandomInt(1, n);
    }
    prev[i] = j;
    return prev;
}

// Generates random instance A of partitioning problem
vector<uint64_t> generateRandomInstance(int n) {
    vector<uint64_t> instance;
    for (int i = 0; i < n; i++) {
        instance.push_back(generateRandomInt(1, pow(10, 12)));
    }
    return instance;
}

// Generates a random sequence solution
vector<int> generateRandomSequenceSoln(int n) {
    vector<int> soln;
    for (int i = 0; i < n; i++) {
        soln.push_back(generateRandomInt(0, 1) == 1 ? 1 : -1);
    }
    return soln;
}

// Generates a random sequence solution MOVE
vector<int> generateRandomSequenceMove(vector<int> prevSeq) {
    int n = prevSeq.size();
    int i = generateRandomInt(0, n - 1);
    int j = generateRandomInt(0, n - 1);
    while (i == j) {
        i = generateRandomInt(0, n - 1);
        j = generateRandomInt(0, n - 1);
    }

    prevSeq[i] = -1 * prevSeq[i];

    // With probability of 0.5, set s_j to -s_j
    int randProb = generateRandomInt(0,1);
    if (randProb == 0) {
        prevSeq[j] = -1 * prevSeq[j];
    }

    return prevSeq;
}

// Given problem instance A and sequence solution S, calculates the residue
uint64_t calcResidueSequence(vector<uint64_t> A, vector<int> seq) {
    int n = seq.size();
    uint64_t set1 = 0;
    uint64_t set2 = 0;
    for (int i = 0; i < n; i++) {
        if (seq[i] == 1) {
            set1 += A[i];
        } else {
            set2 += A[i];
        }
    }
    
    if (set1 > set2) {
        return set1 - set2;
    } else {
        return set2 - set1;
    }
}

// Generates random 64 bit integer between low and high inclusive
uint64_t generateRandomInt(uint64_t low, uint64_t high) {
    std::uniform_int_distribution<uint64_t> distr(low, high);
    return distr(gen);
}

