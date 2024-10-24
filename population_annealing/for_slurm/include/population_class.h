#include <iostream>
#include <fstream>
#include <vector>
#include <bitset>
#include <memory>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
// #include <omp.h>
#include <stdio.h>
#include "/opt/homebrew/Cellar/libomp/18.1.8/include/omp.h" // For parallelizing
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
    double overlap;
    int avg_cluster_size, avg_nowrap_cluster_size;
    int wrap_counter, nowrap_count;


    void reSample(double *Beta, gsl_rng *r, double avg_e, double var_e);
    void calculateEnergies(double *avg_e, double *var_e);
    void calculateFamilies(void);

public:
    // Constructors
    Population(void);
    Population(int nom_pop, gsl_rng *r, double kappa, string mode);
    
    
    void run(string);
    void runSA(string kappastr);
    void runTR(string kappastr);
    void doTwoReplica(double *Beta, fftw_plan p, int num_steps);
    void doTwoRepStep(double padd1, double padd2, fftw_plan p, int num_steps, int m);
    void collectData(double *Beta, double, double);
    void collectDataSA(double *Beta);
    void loadData(string);
    void measureOverlap();

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
};
