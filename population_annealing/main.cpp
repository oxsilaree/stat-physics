#include <iostream>
#include <fstream>
#include <vector>
#include <bitset>
#include <math.h>
#include <memory>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <stdio.h>
//#include <omp.h>
#include "/opt/homebrew/Cellar/libomp/18.1.8/include/omp.h" // For parallelizing
// #include "./lib/libomp/18.1.8/include/omp.h" // For parallelizing in SLURM
#include <time.h>
#include <chrono>
#include <stack>
#include <random>
#include <list>
#include <vector>
#include <stdexcept>

// #include "spin_class.h"
// #include "lattice_class.h"
#include "./include/population_class.h"
#include "./include/functions.h"

using namespace std;

int main(int argc, char** argv)
{
// -------- Catch errors in input arguments

try 
{
    if (argc <= 2) 
    {
        throw std::runtime_error("Kappa/mode not provided. Please type argument in command line.\nTypical values: 0 < kappa < 2\nAll modes: 't', 's', 'p'");
    }
}   catch (const std::runtime_error& e) 
    {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1; 
    }


// -------- Define some constants
    // static int neighbor_table[LEN][LEN][NN_MAX][DIM]; // 6 neighbors (2D ANNNI), 2 coordinates
    gsl_rng *r;
    // int* p;
    // double T;
    int seed = 12345; // We can make this an input later
    srand(time(NULL));
    string prekappastr = argv[1];
    string mode = argv[2];
    double prekappa = stod(argv[1]);
    double kappa = prekappa;                // FOR TESTING
    // double kappa = prekappa/4;          // For N kappa values, we divide by N-1 (FOR BATCH ARRAY JOB)
    string kappastr = to_string(kappa);
    kappastr = kappastr.substr(0,4);    // Precision to 2 decimal places

// -------- Actual code
    initializeRNG(&r, seed);
    auto start = chrono::high_resolution_clock::now(); // for checking time of run
    

    // Population Annealing
    makeNeighborTable();
    cout << "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n";
    cout << "Neighbor table created. \n";
    cout << "Starting simulation...\n";
    cout << "kappa = " << kappa << ", L = " << LEN << ".\n";
    cout << "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n";


    if (mode == "p"){
        Population test_pop(INIT_POP_SIZE, r, kappa, mode);
        cout << "Population Annealing\n";
        cout << "Starting population size = " << INIT_POP_SIZE << ".\n";
        test_pop.run(kappastr);
    } else if (mode == "s") {
        Population test_pop(1, r, kappa, mode);
        cout << "Simulated Annealing\n";
        cout << "No. of Blocks = " << INIT_POP_SIZE << ".\n";
        test_pop.runSA(kappastr);
    } else if (mode == "t") {
        Population test_pop(INIT_POP_SIZE, r, kappa, mode);
        cout << "Two Replica with Annealing\n";
        cout << "Starting population size = " << INIT_POP_SIZE << ".\n";
        test_pop.runTR(kappastr);
    } else {
        cout << "Incorrect mode indicated. Please use the following:\n\
        't' for Two-Replica Cluster with Population Annealing\n\
        's' for Wolff Cluster with Simulated Annealing\n\
        'p' for Wolff Cluster with Population Annealing";
        return 2;
    }
   
    

    auto end = chrono::high_resolution_clock::now(); // For checking duration of program
    chrono::duration<double> elapsed = end - start;

    std::cout << "Time taken: " << elapsed.count() << " seconds." << std::endl;
    cout << "Simulation complete." << endl;
    cout << "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n";
    return 0;
    
}