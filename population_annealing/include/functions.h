#include <iostream>
#include <fstream>
#include <vector>
#include <bitset>
#include <math.h>
#include <memory>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include <time.h>
#include <chrono>
#include <stack>
#include <random>
#include <list>
#include <vector>
#include "parameters.h"

using namespace std;

extern int neighbor_table[LEN][LEN][NN_MAX][DIM];


int* getNeighbor(int, int, int);
void makeNeighborTable();
void initializeRNG(gsl_rng **r, int seed);


