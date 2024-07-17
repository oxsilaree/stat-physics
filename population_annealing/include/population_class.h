#include <iostream>
#include <fstream>
#include <vector>
#include <bitset>
#include <memory>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
// #include <omp.h>
#include <stdio.h>
#include "/usr/local/opt/libomp/include/omp.h" // For parallelizing
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
    int neighbor_table[LEN][LEN][NN_MAX][DIM];
    gsl_rng *r;
    void reSample(double *Beta, gsl_rng *r, double avg_e, double var_e);
    void energy_calcs(double *avg_e, double *var_e);



public:
    Population(void);
    Population(int nom_pop, gsl_rng *r, int nn_table[LEN][LEN][NN_MAX][DIM], double kappa);
    
    
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
    vector<double> beta_values;
    

};
