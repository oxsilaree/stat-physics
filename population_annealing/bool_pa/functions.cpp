#include "functions.h"

void read_NT(int neighbor_table[LENGTH_3][6]) {
	FILE *fp = fopen("neighbors.txt", "r");
	int i = 0;
	while (fscanf(fp, "%d %d %d %d %d %d", &neighbor_table[i][0], &neighbor_table[i][1], &neighbor_table[i][2], &neighbor_table[i][3], &neighbor_table[i][4], &neighbor_table[i][5]) != EOF) {
		i++;
	}
	fclose(fp);
	if (i != LENGTH_3)
		throw "Neighbor table initialization exception.";
}

void read_BT(double bond_table[LENGTH_3][6]) {
	FILE *fp = fopen("bonds.txt", "r");
	int i = 0;
	while (fscanf(fp, "%lf %lf %lf %lf %lf %lf", &bond_table[i][0], &bond_table[i][1], &bond_table[i][2], &bond_table[i][3], &bond_table[i][4], &bond_table[i][5]) != EOF) {
		i++;
	}
	fclose(fp);
	if (i != LENGTH_3)
		throw "Bond table initialization exception.";
}

void initialize_rng(gsl_rng **r, int seed) {
	gsl_rng_env_setup();
	const gsl_rng_type *T;
	T = gsl_rng_mt19937;
	*r = gsl_rng_alloc(T);
	gsl_rng_set(*r, seed);

	for (int i = 0; i < 10000; i++)
		gsl_rng_uniform(*r);
}

int mod(int a, int b) {
	int ret;

	if (b < 0)
		return mod(-a, -b);
	ret = a % b;
	if (ret < 0)
		ret += b;
	return ret;
}

int initialize_NT(int neighbor_table[LENGTH_3][6]) {
	int i, j, k;
	for (i = 0; i < LENGTH; i++) {
		for (j = 0; j < LENGTH; j++) {
			for (k = 0; k < LENGTH; k++) {
				neighbor_table[array_index(k, j, i)][0] = array_index(k, j, mod(i - 1, LENGTH));
				neighbor_table[array_index(k, j, i)][1] = array_index(k, j, mod(i + 1, LENGTH));
				neighbor_table[array_index(k, j, i)][2] = array_index(k, mod(j - 1, LENGTH), i);
				neighbor_table[array_index(k, j, i)][3] = array_index(k, mod(j + 1, LENGTH), i);
				neighbor_table[array_index(k, j, i)][4] = array_index(mod(k - 1, LENGTH), j, i);
				neighbor_table[array_index(k, j, i)][5] = array_index(mod(k + 1, LENGTH), j, i);
			}
		}
	}

	return 0;
}

/* 
For 6 neighbors (3D), there are 2^7=128 local spin configurations. I use the bitset variable because it easily translates between integers and bits.
For each spin location/bonds, this routine iterates over all 128 different possible combinations of spins (bits) and calculates the ratios.
*/
void calc_ratios(double beta, double ratio_table[LENGTH_3][128], double bond_table[LENGTH_3][6]) {
	bitset<7> b;
	int p[2] = {-1, 1};
	for (int s = 0; s < LENGTH_3; s++) {
		for (int i = 0; i < 128; i++) {
			b = i;
			int ss = p[b[0]];
			double sum = 0.0;
			for (int nn = 0; nn < 6; nn++) {
				sum += p[b[nn + 1]] *bond_table[s][nn];
			}
			ratio_table[s][i] = exp(-2 *beta *ss *sum);
		}
	}

}
