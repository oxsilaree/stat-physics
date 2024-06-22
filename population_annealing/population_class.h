#include <iostream>
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
    unique_ptr<Lattice[]> pop_array;
    int nom_pop, max_pop, pop_size ;
    double padd1, padd2;
    double rho_t;
    int neighbor_table[LEN][LEN][nn_max][dim];
    gsl_rng *r;

public:
    Population(void);
    Population(int nom_pop, gsl_rng *r, int nn_table[LEN][LEN][nn_max][dim]);
    void reSample(double T, gsl_rng *r);
    void run(void);

};