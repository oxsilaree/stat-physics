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
    int lattice_number; // Relevant for population annealing
    static spinSite lattice_object[L][L]; // 2 Dimensional now

    // Methods
    void initializeSites();
    void doBurnIn(int neighbor_table[L][L][nn_max][dim]);
    void doBurnInStep(int neighbor_table[L][L][nn_max][dim], double padd1, double padd2);
    void doStep(int neighbor_table[L][L][nn_max][dim], double padd1, double padd2);
    void doSweep(int neighbor_table[L][L][nn_max][dim], double T);
    void doWolffAlgo(int neighbor_table[L][L][nn_max][dim], double T);

    // Update data members
    void updateTotalEnergy(int neighbor_table[L][L][nn_max][dim]);
    void updateTotalMag();

    // Get data members
    int getTotalEnergy();
    int getTotalMag();
};