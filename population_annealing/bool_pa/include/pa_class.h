#include <iostream>
#include <vector>
#include <bitset>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include <omp.h>
#include <memory>
#include "params.h"
#include "functions.h"
#include "ising_class.h"

using namespace std;

class PA_simulation {
	private:
		unique_ptr<spin_conf[]> pop_array;
		int nom_pop;
		int max_pop;
		int pop_size;
		gsl_rng *r;
		int unique_families;
		double rho_t;
		double gs_e;
		spin_conf gs;
		int neighbor_table[LENGTH_3][6];
		double bond_table[LENGTH_3][6];
		double ratio_table[LENGTH_3][128];
		void resample(double *beta, double avg_e, double var_e, gsl_rng *r);
		void energy_calcs(double *avg_e, double *var_e);
		void family_calcs(void);
	public:
		PA_simulation(void);
		PA_simulation(int nom_pop, gsl_rng *r, int neighbor_table[LENGTH_3][6], double bond_table[LENGTH_3][6]);
		void run(void);
		spin_conf get_gs(void);
		double get_gs_e(void);
		double get_rho_t(void);
		int num_gs_families(void);
};
