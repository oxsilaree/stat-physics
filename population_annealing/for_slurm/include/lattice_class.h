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
    int mag, abs_mag, avg_cluster_size, avg_nowrap_cluster_size;
    int wrap_counter, nowrap_counter; // Relevant for data collection
    double energy, spec_heat, suscep;  //---^
    double kappa;
    vector<vector<spinSite> > lattice_object;
    // mutex lattice_mutex;
    

public:
    // Public member attribute (so we can change/check it easily)
    // spinSite lattice_object[LEN][LEN]; // 2 Dimensional now
    

    int lattice_number; // Relevant for population annealing
    
    // Constructor Declaration
    Lattice(double kappa);

    // Methods
    void initializeSites(int neighbor_table[LEN][LEN][NN_MAX][DIM], double *Beta);
    void doBurnIn(int neighbor_table[LEN][LEN][NN_MAX][DIM], double Beta);
    void doBurnInStep(int neighbor_table[LEN][LEN][NN_MAX][DIM], double *padd1, double *padd2);
    void doStep(int neighbor_table[LEN][LEN][NN_MAX][DIM], double *padd1, double *padd2);
    void doSweep(int neighbor_table[LEN][LEN][NN_MAX][DIM], double *Beta);
    void doWolffAlgo(int neighbor_table[LEN][LEN][NN_MAX][DIM], double *Beta);

    // Update data members
    void updateTotalEnergy(int neighbor_table[LEN][LEN][NN_MAX][DIM]);
    void updateTotalMag();

    // Get data members
    int getTotalEnergy();
    int getTotalMag();
    double getAvgClusterSize();
    double getAvgNowrapClusterSize();
    int getNoWrapCount();
    spinSite* getSpinSite(int row, int col);

    // spinSite(&getLattice())[LEN][LEN] { return lattice_object; };
};

inline double Lattice::getAvgClusterSize()
{
    return avg_cluster_size;
}

inline double Lattice::getAvgNowrapClusterSize()
{
    return avg_nowrap_cluster_size;
}

inline int Lattice::getNoWrapCount()
{
    return nowrap_counter;
}


inline int Lattice::getTotalEnergy()
{
    return energy;
}

inline int Lattice::getTotalMag()
{
    return mag;
}

inline spinSite* Lattice::getSpinSite(int row, int col)
{
    // cout <<"---row = " << row << ", col = " << col << "---\n";
    return &lattice_object[(int)row][(int)col];
}