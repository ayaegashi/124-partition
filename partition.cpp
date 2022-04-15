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

using namespace std;

int kk(vector<int> seq);

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

    // Construct vector from inputfile
    vector<int> sequence;
    ifstream file (inputfile, ios::in);
    if (file.fail()) {
        printf("File could not be opened.\n");
        return 0;
    }
    int n;
    while (file >> n){
        sequence.push_back(n);
    }

    file.close();

    for (size_t i = 0; i < sequence.size(); i ++) {
        printf("%i ", sequence[i]);
    }

    if (flag == 0 && alg == 0){
        int difference = kk(sequence);
        printf("\n%i\n", difference);
    }

    return 0;
};

// IGNORE FOR NOW DOESNT WORK!
// Karmarkar-Karp Algorithm --> returns just difference
int kk(vector<int> seq){
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