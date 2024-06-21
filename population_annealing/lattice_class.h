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
#include "parameters.h"
#include "functions.h"
#include "spin_class.h"

using namespace std;

class Lattice
{
private:
    int energy, mag, abs_mag, avg_cluster_size; // Relevant for data collection
    double spec_heat, suscep; //--^
    
    

public:
    // Public member attribute (so we can change/check it easily)
    // spinSite lattice_object[LEN][LEN]; // 2 Dimensional now
    unique_ptr<unique_ptr<spinSite[]>[]> lattice_object;
    int lattice_number; // Relevant for population annealing
    

    // Constructor

    Lattice()
    {
        lattice_object = std::make_unique<std::unique_ptr<spinSite[]>[]>(LEN);
        for (int i = 0; i < LEN; ++i) {
            lattice_object[i] = std::make_unique<spinSite[]>(LEN);
        }
    }

    // Methods
    void initializeSites();
    void doBurnIn(int neighbor_table[LEN][LEN][nn_max][dim], double T);
    void doBurnInStep(int neighbor_table[LEN][LEN][nn_max][dim], double padd1, double padd2);
    void doStep(int neighbor_table[LEN][LEN][nn_max][dim], double padd1, double padd2);
    void doSweep(int neighbor_table[LEN][LEN][nn_max][dim], double T);
    void doWolffAlgo(int neighbor_table[LEN][LEN][nn_max][dim], double T, int num_sweeps);

    // Update data members
    void updateTotalEnergy(int neighbor_table[LEN][LEN][nn_max][dim]);
    void updateTotalMag();

    // Get data members
    int getTotalEnergy();
    int getTotalMag();

    // spinSite(&getLattice())[LEN][LEN] { return lattice_object; };
};