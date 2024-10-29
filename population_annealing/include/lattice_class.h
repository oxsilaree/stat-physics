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
#include <fftw3.h>
#include <deque>
#include <set>

using namespace std;

class Lattice
{
private:
    int mag, abs_mag;
    double energy, spec_heat, suscep;  //---^
    double kappa;
    double dom_freq, dom_amplitude; // Dominant frequency and corresponding amplitude
    vector<vector<spinSite> > lattice_object;
    int family, new_family; // For checking ancestors of replicas (and for pairing)
    // std::deque<int> recent_families;

    // Wrapping and percolation observables
    int wrap_counter, no_wrap_counter;
    int z_wrap_counter, x_wrap_counter, xz_wrap_counter;
    double avg_cluster_size, avg_nowrap_cluster_size;
    double avg_zwrap_cluster_size, avg_xwrap_cluster_size, avg_xzwrap_cluster_size;

   
    

public:
    // Public member attribute (so we can change/check it easily)
    // spinSite lattice_object[LEN][LEN]; // 2 Dimensional now
    

    int lattice_number; // Relevant for population annealing
    
    
    Lattice(double kappa, int family); // Constructor
    //  Lattice(const Lattice&); // Copy constructor

    // Methods
    void initializeSites(double *Beta);
    void doBurnIn(double Beta);
    void doBurnInStep(double *padd1, double *padd2);
    void doStep(double *padd1, double *padd2);
    void doSweep(double *Beta);
    void doWolffAlgo(double *Beta, fftw_plan p, int num_steps);
    void doFFT(fftw_plan p);
    

    // Update data members
    void updateTotalEnergy();
    void updateTotalMag();

    // Get data members
    int getTotalEnergy();
    int getTotalMag();
    double getAvgClusterSize();
    double getAvgNoWrapClusterSize();
    double getAvgXWrapClusterSize();
    double getAvgZWrapClusterSize();
    double getAvgXZWrapClusterSize();
    int getWrapCount();
    int getNoWrapCount();
    int getZWrapCount();
    int getXWrapCount();
    int getXZWrapCount();
    double getDomFreq();
    double getDomAmplitude();
    spinSite* getSpinSite(int row, int col);
    int getFamily();
    int getNewFamily();
    void setNewFamily(int nf);
    deque<int> getRecentFamilies();
    void updateRecentFamilies(int nf);
    void printLattice();

};

inline double Lattice::getAvgClusterSize()
{
    return avg_cluster_size;
}

inline double Lattice::getAvgNoWrapClusterSize()
{
    return avg_nowrap_cluster_size;
}

inline double Lattice::getAvgXWrapClusterSize()
{
    return avg_xwrap_cluster_size;
}

inline double Lattice::getAvgZWrapClusterSize()
{
    return avg_zwrap_cluster_size;
}

inline double Lattice::getAvgXZWrapClusterSize()
{
    return avg_xzwrap_cluster_size;
}


inline int Lattice::getWrapCount()
{
    return wrap_counter;
}

inline int Lattice::getNoWrapCount()
{
    return no_wrap_counter;
}

inline int Lattice::getZWrapCount()
{
    return z_wrap_counter;
}

inline int Lattice::getXWrapCount()
{
    return x_wrap_counter;
}

inline int Lattice::getXZWrapCount()
{
    return xz_wrap_counter;
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

inline double Lattice::getDomFreq()
{
    return dom_freq;
}

inline double Lattice::getDomAmplitude()
{
    return dom_amplitude;
}

inline int Lattice::getFamily()
{
    return family;
}

inline int Lattice::getNewFamily()
{
    return new_family;
}

inline void Lattice::setNewFamily(int nf)
{
    new_family = nf;
}

/*
inline deque<int> Lattice::getRecentFamilies()
{
    return recent_families;
}*/