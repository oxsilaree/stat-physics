#include <iostream>
#include <vector>
#include <bitset>
#include <math.h>
#include <memory>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <chrono>
#include <stack>
#include <random>
#include <list>
#include <vector>
// #include "spin_class.h"
#include "lattice_class.h"
// #include "population_class.h"
#include "functions.h"

using namespace std;

int main()
{
// -------- Constants
    static int neighbor_table[L][L][nn_max][dim]; // 6 neighbors (2D ANNNI), 2 coordinates
    int* p;

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


// Make neighbor table (this works CAA 17/6/24)
    for (int i = 0; i < L; i++)
    {
        for (int j = 0; j < L; j++)
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
    test_lattice.doWolffAlgo(neighbor_table, 10);
    for (int i = 0; i < L; i++)
    {
        for (int j = 0; j < L; j++)
        {
            cout << test_lattice.lattice_object[i][j].getSpin() + 3 << ", ";
        }
        cout << endl;
    }
    return 0;
}