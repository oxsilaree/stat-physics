#include <iostream>
#include <fstream>
#include <vector>
#include <bitset>
#include <memory>
#include <gsl/gsl_rng.h> // GNU Standard Library
#include <gsl/gsl_randist.h>
#include <stdio.h>
#include <omp.h> // OpenMP for parallelizing
#include <math.h>
#include <time.h>
#include <chrono>
#include <stack>
#include <random>
#include <list>
#include <vector>
#include "parameters.h"
#include "functions.h"
#include "lattice_class.h"

using namespace std;

class Population
{
private:
    // unique_ptr<Lattice[]> pop_array;
    vector<Lattice> pop_array;
    int nom_pop, max_pop, pop_size;
    double padd1, padd2;
    double rho_t;
    double kappa;
    gsl_rng *r;
    void reSample(double *Beta, gsl_rng *r, double avg_e, double var_e);
    void energy_calcs(double *avg_e, double *var_e);



public:
    Population(void);
    Population(int nom_pop, gsl_rng *r, double kappa);
    
    
    void run(string);
    void collectData(double *Beta, double, double);
    void loadData(string);

    // double energy_data[(int)T_ITER];
    // double spec_heat_data[(int)T_ITER];
    // double magnetization_data[(int)T_ITER];
    // double susceptibility_data[(int)T_ITER];

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
};
