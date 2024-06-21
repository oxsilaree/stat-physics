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
#include "lattice_class.h"
// #include "population_class.h"
#include "functions.h"

using namespace std;

int main(int argc, char** argv)
{
// -------- Constants
    static int neighbor_table[LEN][LEN][nn_max][dim]; // 6 neighbors (2D ANNNI), 2 coordinates
    int* p;
    double T;
    srandom(time(NULL)); // Use consistent random seed to make sure behavior is consistent
    // double kappa = stod(argv[1]);

// -------- Lists for data
    list<double> E;     //-| v
    list<double> M;     //-| v
    list<double> Mabs;  //-| For typical data
    list<double> C;     //-| ^
    list<double> X;     //-| ^
    list<double> Mz; // For order parameter (modulated waveform)
    list<int> cluster_sizes; // For percolation

// -------- Catch errors in input arguments



// -------- Actual code

    cout << "kappa = " << kappa << endl;
// Make neighbor table (this works CAA 17/6/24)
    for (int i = 0; i < LEN; i++)
    {
        for (int j = 0; j < LEN; j++)
        {
            for (int pos = 0; pos < nn_max; pos++)
            {
                p = getNeighbor(i, j, pos); 
                neighbor_table[i][j][pos][0] = p[0];
                neighbor_table[i][j][pos][1] = p[1];
                // cout << neighbor_table[i][j][pos][0] << ", " << neighbor_table[i][j][pos][1] << endl;  
            }
        }
    }

    /* TEST STUFF OUT */

    Lattice test_lattice;
    test_lattice.initializeSites();
    for (int i = 0; i < LEN; i++)
    {
        for (int j = 0; j < LEN; j++)
        {
            cout << test_lattice.lattice_object[i][j].getSpin() + 2 << ", ";
        }
        cout << "\n";
    }
    cout << ".\n.\n.\n";
    
    // ofstream test_data;
    // test_data.open("/Users/shanekeiser/Documents/Summer 2024/Research/PopulationAnnealing/data/test_data");
    cout << "Hello world!" << endl;
    T = T_init;
    for (int l = 0; l < T_iter; l++) {
        test_lattice.doWolffAlgo(neighbor_table, T, SWEEPS);
        cout << "Finished run for T = " << T << "." << endl;
        T -= double(T_init/T_iter); // "Cooling" the system
        T = floor((100.*T)+.5)/100;
    }
    
    for (int i = 0; i < LEN; i++)
    {
        for (int j = 0; j < LEN; j++)
        {
            cout << test_lattice.lattice_object[i][j].getSpin() + 2 << ", ";
            // test_data << test_lattice.lattice_object[i][j].getSpin() + 2 << ", ";
        }
        cout << "\n";
        // test_data << "\n";
    }
    // test_data.close();
    cout << "Simulation complete." << endl;
    return 0;
}