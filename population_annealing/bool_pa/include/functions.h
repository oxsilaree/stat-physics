#include <iostream>
#include <vector>
#include <bitset>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include <omp.h>
#include <memory>
#include "params.h"

using namespace std;

void read_NT(int neighbor_table[LENGTH_3][6]);
void read_BT(double bond_table[LENGTH_3][6]);
void initialize_rng(gsl_rng **r, int seed);
int mod(int a, int b);
int initialize_NT(int neighbor_table[LENGTH_3][6]);
void calc_ratios(double beta, double ratio_table[LENGTH_3][128], double bond_table[LENGTH_3][6]);
