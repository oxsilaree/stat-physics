#include "pa_class.h"

/*
The PA_simulation initializer requires inputs and otherwise returns an error.
*/
PA_simulation::PA_simulation(void) {
	throw invalid_argument("Invalid PA simulation initialization parameters.");
}

/*
The pop_array is a unique pointer that needs to point to a newly allocated array of spin configurations. This is so the
simulation can have a variable population size. Unique pointers will automatically deallocate when they are destroyed, 
so this makes "trash collection" automatic.
*/
PA_simulation::PA_simulation(int nom_pop, gsl_rng *r, int neighbor_table[LENGTH_3][6], double bond_table[LENGTH_3][6]) {
	PA_simulation::nom_pop = nom_pop;
	PA_simulation::max_pop = nom_pop + sqrt(nom_pop *10);
	PA_simulation::pop_size = nom_pop;

	PA_simulation::pop_array = unique_ptr<spin_conf[]>(new spin_conf[PA_simulation::max_pop]);

	PA_simulation::r = r;	

	for (int s = 0; s < LENGTH_3; s++) {
		for (int n = 0; n < 6; n++) {
			PA_simulation::neighbor_table[s][n] = neighbor_table[s][n];
			PA_simulation::bond_table[s][n] = bond_table[s][n];
		}
	}

	for (int m = 0; m < nom_pop; m++) {
		pop_array[m].init(r, m, neighbor_table, bond_table);
		pop_array[m].calc_energy(neighbor_table, bond_table);
	}

	double dummy_1, dummy_2;
	energy_calcs(&dummy_1, &dummy_2);
}

void PA_simulation::family_calcs(void) {
	int i;
	int family_size_dist[nom_pop];
	int family_size[nom_pop];

	for (i = 0; i < nom_pop; i++) {
		family_size[i] = 0;
		family_size_dist[i] = 0;
	}

//	Count the number of ancestors for each family
	for (i = 0; i < pop_size; i++)
		family_size[pop_array[i].get_family()]++;

//	Count the distribution of family sizes
	for (i = 0; i < nom_pop; i++)
		family_size_dist[pop_array[i].get_family()]++;

	int counter = 0;
	rho_t = 0;

	for (i = 0; i < nom_pop; i++) {
		if (family_size[i] > 0) {
			counter++;
			rho_t += (double) family_size[i] *family_size[i] /pop_size /pop_size;
		}
	}

	unique_families = counter;
	rho_t *= nom_pop;
}

int PA_simulation::num_gs_families(void) {
	int i;
	int family_size[nom_pop];

	for (i = 0; i < nom_pop; i++) {
		family_size[i] = 0;
	}

//	Count the number of ancestors for each family
	for (i = 0; i < pop_size; i++)
		if (pop_array[i].get_energy() == gs_e)
			family_size[pop_array[i].get_family()]++;

	int counter = 0;
	for (i = 0; i < nom_pop; i++) {
		if (family_size[i] > 0)
			counter++;
	}

	return counter;
}

void PA_simulation::resample(double *beta, double avg_e, double var_e, gsl_rng *r) {
	int m, new_pop_size = 0;
	double config_weight[pop_size];
	double Q = 0;
	int num_replicas[pop_size];

	double delta_beta = CULLING_FRAC *sqrt(2 *PI /var_e);

	if (*beta + delta_beta < 5)
		*beta += delta_beta;
	else
		*beta = 5.0;

//	Calculate configuration weight
	for (m = 0; m < pop_size; m++) {
		double temp_e = pop_array[m].get_energy();
		double temp_w = -(delta_beta) *(temp_e - avg_e);
		config_weight[m] = exp(temp_w);
		Q += config_weight[m];
	}

	Q /= (double) pop_size;

//	Calculate the number of replicas (probabilistically)
	for (m = 0; m < pop_size; m++) {
		double tau = (nom_pop / ((double) pop_size)) *(config_weight[m] /Q);
		int floor = (int) tau;
		int ceiling = floor + 1;

		if (gsl_rng_uniform(r) < ceiling - tau)
			num_replicas[m] = floor;
		else
			num_replicas[m] = ceiling;

		new_pop_size += num_replicas[m];
	}

	if (new_pop_size > max_pop)
		throw "Maximum population size exceeded.";
/*
	Fill empty gaps in population array so un-erased members are contiguous
	(e.g. this string corresponds to the population and the numbers correspond to the number of
	copies of each replica, then this routine changes 12002004013 into 12312400000)

	Instead of this complicated routine, we could use vectors with their delete function, but
	vectors can't be parallelized using omp and deleting is slow.
*/
	int copy_to = 0, copy_from = pop_size - 1;
	while (copy_to < copy_from) {
		while (num_replicas[copy_to] > 0)
			copy_to++;
		while (num_replicas[copy_from] <= 0)
			copy_from--;
		if (copy_to < copy_from) {
			pop_array[copy_to] = pop_array[copy_from];
			num_replicas[copy_to] = num_replicas[copy_from];
			num_replicas[copy_from] = 0;
		}
	}

	pop_size = new_pop_size;
/*
	Make copies of config's from the beginning to the end of contiguous data block
	(12312400000 and ABCDEFGHIJK -> ABCDEFBCCEFFF)
*/
	copy_to = 0;
	while (num_replicas[copy_to] > 0)
		copy_to++;
	int copy_end = copy_to;
	copy_from = 0;
	while (copy_from < copy_end) {
		for (m = 0; m < num_replicas[copy_from] - 1; m++) {
			pop_array[copy_to] = pop_array[copy_from];
			copy_to++;
		}
		copy_from++;
	}
}

void PA_simulation::run(void) {
	gsl_rng *r_thread[NUM_THREADS];
	double avg_e = 0.0, var_e = 0.0;
	FILE *fp = fopen("run.dat", "w");

	for (int i = 0; i < NUM_THREADS; i++) {
		int tmp_rnd = gsl_rng_uniform_int(r, 10000000) + 1000000;
		initialize_rng(&r_thread[i], tmp_rnd);
	}
	double beta = 0.0;
	int num_sweeps;

	while (beta < MAX_BETA) {
		fprintf(fp, "%5.5f\n", beta);
		fflush(fp);
		if (beta < 0.5)
			num_sweeps = 3;
		else
			num_sweeps = 21;
		energy_calcs(&avg_e, &var_e);
		resample(&beta, avg_e, var_e, r);
		calc_ratios(beta, ratio_table, bond_table);
		if (beta < 0.1) {
			#pragma omp parallel for shared(pop_array, neighbor_table, bond_table, beta, r_thread, ratio_table)
			for (int m = 0; m < pop_size; m++) {
				int thread = omp_get_thread_num();
				pop_array[m].high_T_Met_sweep_set(r_thread[thread], neighbor_table, bond_table, ratio_table, num_sweeps);
				pop_array[m].calc_energy(neighbor_table, bond_table);
			}
		}
		else {
			#pragma omp parallel for shared(pop_array, neighbor_table, bond_table, beta, r_thread, ratio_table, num_sweeps)
			for (int m = 0; m < pop_size; m++) {
				int thread = omp_get_thread_num();
				pop_array[m].Met_sweep_set(r_thread[thread], neighbor_table, bond_table, ratio_table, num_sweeps);
				pop_array[m].calc_energy(neighbor_table, bond_table);
			}
		}
	}
	fclose(fp);
	energy_calcs(&avg_e, &var_e);
	family_calcs();
}

void PA_simulation::energy_calcs(double *avg_e, double *var_e) {
	int min_index = 0;
	double temp_min = 0;
	double avg_e2 = 0.0;
	double temp_avg_e = 0.0;
	for (int m = 0; m < pop_size; m++) {
		double temp_e = pop_array[m].get_energy();
		temp_avg_e += temp_e;
		avg_e2 += temp_e *temp_e;
		if (temp_e < temp_min) {
			temp_min = temp_e;
			min_index = m;
		}
	}
	avg_e2 /= pop_size;
	temp_avg_e /= pop_size;
	*avg_e = temp_avg_e;
	*var_e = (avg_e2 - temp_avg_e *temp_avg_e);
	gs = pop_array[min_index];
	gs_e = gs.get_energy();
}

spin_conf PA_simulation::get_gs(void) {
	return gs;
}

double PA_simulation::get_gs_e(void) {
	return gs_e;
}

double PA_simulation::get_rho_t(void) {
	return rho_t;
}
