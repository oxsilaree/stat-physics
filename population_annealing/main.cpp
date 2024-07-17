#include <iostream>
#include <fstream>
#include <vector>
#include <bitset>
#include <math.h>
#include <memory>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
// #include <omp.h>
#include "/usr/local/opt/libomp/include/omp.h" // For parallelizing
#include <time.h>
#include <chrono>
#include <stack>
#include <random>
#include <list>
#include <vector>
#include <stdexcept>

// #include "spin_class.h"
// #include "lattice_class.h"
#include "population_class.h"
#include "functions.h"

using namespace std;

int main(int argc, char** argv)
{
// -------- Catch errors in input arguments

try 
{
    if (argc <= 1) 
    {
        throw std::runtime_error("Kappa not provided. Please type argument in command line.\n(Typical values: 0 < kappa < 2)\n");
    }
}   catch (const std::runtime_error& e) 
    {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1; 
    }



// -------- Define some constants
    static int neighbor_table[LEN][LEN][NN_MAX][DIM]; // 6 neighbors (2D ANNNI), 2 coordinates
    gsl_rng *r;
    int* p;
    // double T;
    int seed = 1; // We can make this an input later
    srand(time(NULL));
    string kappastr = argv[1];
    double kappa = stod(argv[1]);


// -------- Lists for data
    list<double> E;     //-| v
    list<double> M;     //-| v
    list<double> Mabs;  //-| For typical data
    list<double> C;     //-| ^
    list<double> X;     //-| ^
    list<double> Mz; // For order parameter (modulated waveform)
    list<int> cluster_sizes; // For percolation

// -------- Actual code
    initializeRNG(&r, seed);
    auto start = chrono::high_resolution_clock::now(); // for checking time of run
    
    for (int i = 0; i < LEN; i++) // Make neighbor table
    {
        for (int j = 0; j < LEN; j++)
        {
            for (int pos = 0; pos < NN_MAX; pos++)
            {
                p = getNeighbor(i, j, pos); 
                neighbor_table[i][j][pos][0] = p[0];
                neighbor_table[i][j][pos][1] = p[1];
                // cout << neighbor_table[i][j][pos][0] << ", " << neighbor_table[i][j][pos][1] << endl;  
            }
        }
    }
    
    // Population Annealing

    cout << "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n";
    cout << "Neighbor table created. \n";
    cout << "Starting simulation...\n";
    cout << "kappa = " << kappa << ".\n";
    cout << "Starting population size = " << INIT_POP_SIZE << ".\n";
    cout << "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n";

    omp_set_num_threads(NUM_THREADS);
    Population test_pop(INIT_POP_SIZE, r, neighbor_table, kappa);
    test_pop.run(kappastr);
    

    auto end = chrono::high_resolution_clock::now(); // For checking duration of program
    chrono::duration<double> elapsed = end - start;

    std::cout << "Time taken: " << elapsed.count() << " seconds." << std::endl;
    cout << "Simulation complete." << endl;
    return 0;
    
}