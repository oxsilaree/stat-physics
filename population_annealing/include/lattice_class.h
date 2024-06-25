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
#include <mutex>

using namespace std;

class Lattice
{
private:
    int mag, abs_mag, avg_cluster_size; // Relevant for data collection
    double energy, spec_heat, suscep; //--^
    vector<vector<spinSite> > lattice_object;
    // mutex lattice_mutex;
    

public:
    // Public member attribute (so we can change/check it easily)
    // spinSite lattice_object[LEN][LEN]; // 2 Dimensional now
    

    int lattice_number; // Relevant for population annealing
    
    // Constructor Declaration
    Lattice();

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
    spinSite* getSpinSite(int row, int col);

    // spinSite(&getLattice())[LEN][LEN] { return lattice_object; };
};