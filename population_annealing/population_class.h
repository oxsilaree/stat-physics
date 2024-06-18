#include <iostream>
#include <vector>
#include <bitset>
#include <math.h>
#include <memory>
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
public:
    unique_ptr<Lattice[]> pop_array;
    int nom_pop, max_pop, pop_size ;
    double padd1, padd2;
    double rho_t;
    int neighbor_table[L][L][nn_max][dim];

private:
    Population(void);
    Population(int nom_pop, int neighbor_table[L][L][nn_max][dim]);
};