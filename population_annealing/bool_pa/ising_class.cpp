#include "ising_class.h"


void spin_conf::calc_energy(int neighbor_table[LENGTH_3][6], double bond_table[LENGTH_3][6]) {
	energy = 0.0;
	int p[2] = {-1, 1};
	for (int s = 0; s < LENGTH_3; s++) {
		double temp_sum = 0;
		for (int n = 0; n < 6; n++) {
			int neighbor = neighbor_table[s][n];
			temp_sum += p[spin[neighbor]] *bond_table[s][n];
		}
		energy += -temp_sum *p[spin[s]];
	}
	energy /= 2.0;
}

void spin_conf::print(void) {
	FILE *fp = fopen("conf.dat", "w");
	for (int i = 0; i < LENGTH_3; i++)
		fprintf(fp, "%d\n", (int) spin[i]);
}
void spin_conf::init(gsl_rng *r, int family, int neighbor_table[LENGTH_3][6], double bond_table[LENGTH_3][6]) {
	for (int s = 0; s < LENGTH_3; s++) {
		if (gsl_rng_uniform(r) < 0.5)
			spin[s] = 0;
		else
			spin[s] = 1;
	}
	spin_conf::family = family;
	calc_energy(neighbor_table, bond_table);
}

void spin_conf::ferro_init(gsl_rng *r, int family, int neighbor_table[LENGTH_3][6], double bond_table[LENGTH_3][6]) {
	for (int s = 0; s < LENGTH_3; s++) {
		spin[s] = 1;
	}
	spin_conf::family = family;
	calc_energy(neighbor_table, bond_table);
}

/*
There are two metropolis algorithm routines, the high-T routine randomly chooses a spin so that we avoid the scenario of the whole lattice flipping.
The index variable calculates an integer between 0 and 127 that determines the encoding for the ratio table, everything else is self-explanatory.
*/
double spin_conf::high_T_Met_sweep_set(gsl_rng *r, int neighbor_table[LENGTH_3][6], double bond_table[LENGTH_3][6], double ratio_table[LENGTH_3][128], int sweeps) {
	double acc_ratio = 0.0;
	for (int swp = 0; swp < sweeps; swp++) {
		for (int i = 0; i < LENGTH_3; i++) {
			int s = gsl_rng_uniform_int(r, LENGTH_3);
			int index = spin[s] + 2 *spin[neighbor_table[s][0]] + 4 *spin[neighbor_table[s][1]] + 8 *spin[neighbor_table[s][2]] + 16 *spin[neighbor_table[s][3]] + 32 *spin[neighbor_table[s][4]] + 64 *spin[neighbor_table[s][5]];
			double weight = ratio_table[s][index];
			if (weight >= 1 || gsl_rng_uniform(r) < weight) {
				spin[s] = !spin[s];
				acc_ratio++;
			}
		}
	}
	acc_ratio /= ((double) LENGTH_3 *sweeps);
	return acc_ratio;
}

double spin_conf::Met_sweep_set(gsl_rng *r, int neighbor_table[LENGTH_3][6], double bond_table[LENGTH_3][6], double ratio_table[LENGTH_3][128], int sweeps) {
	double acc_ratio = 0.0;
	for (int swp = 0; swp < sweeps; swp++) {
		for (int s = 0; s < LENGTH_3; s++) {
			int index = spin[s] + 2 *spin[neighbor_table[s][0]] + 4 *spin[neighbor_table[s][1]] + 8 *spin[neighbor_table[s][2]] + 16 *spin[neighbor_table[s][3]] + 32 *spin[neighbor_table[s][4]] + 64 *spin[neighbor_table[s][5]];
			double weight = ratio_table[s][index];
			if (weight >= 1 || gsl_rng_uniform(r) < weight) {
				spin[s] = !spin[s];
				acc_ratio++;
			}
		}
	}
	acc_ratio /= ((double) LENGTH_3 *sweeps);
	return acc_ratio;
}

int spin_conf::get_family(void) {
	return family;
}

double spin_conf::get_energy(void) {
	return energy;
}
