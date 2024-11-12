#include <iostream>
#include <fstream>
#include <vector>
#include <bitset>
#include <memory>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
// #include <omp.h>
#include <stdio.h>
#include <algorithm>
#include "/opt/homebrew/Cellar/libomp/18.1.8/include/omp.h" // For parallelizing
// include "../lib/libomp/18.1.8/include/omp.h" // For parallelizing in SLURM
#include <math.h>
#include <time.h>
#include <chrono>
#include <stack>
#include <random>
#include <utility>
#include <list>
#include "parameters.h"
#include "functions.h"
#include "lattice_class.h"
#include <deque>
#include <set>

#include <fftw3.h>

using namespace std;

class Population
{
private:
    // unique_ptr<Lattice[]> pop_array;
    vector<Lattice> pop_array;
    int nom_pop, max_pop, pop_size;
    double padd1, padd2;
    double kappa;
    string mode;
    gsl_rng *r;
    int unique_families;
    double rho_t;
    double gs_e;
    int num_steps;
    double overlap, var_overlap, abs_overlap;

    // Wrapping and percolation observables
    int wrap_counter;
    int z_wrap_counter, x_wrap_counter, xz_wrap_counter, no_wrap_counter;
    double avg_cluster_size;
    double avg_zwrap_cluster_size, avg_xwrap_cluster_size, avg_xzwrap_cluster_size, avg_nowrap_cluster_size;

    // Other observables
    double free_energy;
    double smoothed_var_e;


    void reSample(double *Beta, gsl_rng *r, double avg_e, double var_e);
    void calculateEnergies(double *avg_e, double *var_e);
    void calculateFamilies(void);
    void makeHistograms(string kappastr);

public:
    // Constructors
    Population(void);
    Population(int nom_pop, gsl_rng *r, double kappa, string mode);
    
    
    void run(string);
    void runSA(string kappastr);
    void runTR(string kappastr);
    void doTwoReplica(double padd1, double padd2, int num_steps, gsl_rng *r, int index1, int index2);
    void doTwoRepStep(double padd1, double padd2, gsl_rng *r, int index1, int index2);
    void collectData(double *Beta, double, double);
    void getStructureFactorIntensity(string kappastr, int);
    void collectDataSA(double *Beta);
    void loadData(string);
    void measureOverlap();
    void getOverlapDistribution(vector<int>, int, string);
    // bool haveSharedFamily(Lattice* lattice1, Lattice* lattice2);

	int countFamilies(void);
    

    vector<double> energy_data;
    vector<double> spec_heat_data;
    vector<double> magnetization_data;
    vector<double> susceptibility_data;
    vector<double> clustersize_data; // Only non-wrapping clusters
    vector<double> nowrapclustersize_data;
    vector<double> beta_values;
    vector<double> energy_sq_data;
    vector<double> magnetization_sq_data;
    vector<double> magnetization_abs_data;
    vector<double> wrapping_data;
    vector<double> fft_freq_data;
    vector<double> fft_amp_data;
    vector<double> rho_t_data;
    vector<double> unique_families_data;
    vector<double> overlap_data;
    vector<double> abs_overlap_data;
    vector<double> var_overlap_data;
    vector<double> free_energy_data;

    vector<double> z_wrapping_data;
    vector<double> x_wrapping_data;
    vector<double> xz_wrapping_data;
    vector<double> z_clustersize_data;
    vector<double> x_clustersize_data;
    vector<double> xz_clustersize_data;

};
