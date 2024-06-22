#include "population_class.h"

using namespace std;

/*
    This population annealing routine has been adapted with much gratitude 
    from the work of C. Amey, a previous student of my professor.

    I have made some slight changes, mostly in naming conventions.
*/

// Initializers
Population::Population(void) 
{
	throw invalid_argument("Invalid PA simulation initialization parameters.");
}

Population::Population(int nom_pop, gsl_rng *r, int nn_table[LEN][LEN][nn_max][dim])
{
    // Initialize population
    Population::nom_pop = nom_pop;
    Population::max_pop = int(sqrt(nom_pop) * 10);
    Population::pop_size = nom_pop;
    Population::pop_array = unique_ptr<Lattice[]>(new Lattice[Population::max_pop]); // Not sure how this works but this makes
                                                                                     // the population
    for (int i = 0; i < LEN; ++i) {
        for (int j = 0; j < LEN; ++j) {
            for (int pos = 0; pos < nn_max; ++pos) {
                for (int d = 0; d < dim; ++d) {
                    Population::neighbor_table[i][j][pos][d] = nn_table[i][j][pos][d];
                }
            }
        }
    }

    Population::r = r;	                    
}

void Population::reSample(double T, gsl_rng *r)
{
    int j, new_pop_size,    floor, ceiling; //= 0;
    double energy_j, weight_j, tau_j, expected_copies_j; // used for reweighting
    double Beta = 1/T;
    double Beta_prime = 1/(T - (T_init-T_final)/T_iter);
    double d_Beta = Beta_prime - Beta; // Get the difference in previous and new inverse temperature
	double config_weight[pop_size];
	double Q = 0;
	int num_replicas[pop_size];

    for (j = 0; j < pop_size; j++) // Get Q(beta, beta') for each lattice in population
    {
        energy_j = pop_array[j].getTotalEnergy();
        weight_j = -d_Beta * energy_j;
        config_weight[j] = exp(weight_j); // Storing these values is good to get tau later on
        Q += config_weight[j];
    }
    Q /= (double)pop_size;

    for (j = 0; j < pop_size; j++) // Get the expected new number of replicas per currently existing lattice,
    {                              // and new population size
        tau_j = config_weight[j]/Q;
        expected_copies_j = (nom_pop/(double)pop_size) * tau_j;
        floor = (int)expected_copies_j;
        ceiling = floor + 1;
        if (gsl_rng_uniform(r) < ceiling - tau_j)
            num_replicas[j] = floor;
        else
            num_replicas[j] = ceiling;
        new_pop_size += num_replicas[j];
    }
    /*
    if (new_pop_size > max_pop)
		throw "Maximum population size exceeded."; // Must 'catch' errors like this later on
    */

    /*
        Fill empty gaps in population array so un-erased members are contiguous
        (e.g. this string corresponds to the population and the numbers correspond to the number of
        copies of each replica, then this routine changes 12002004013 into 12312400000)
    */

	int copy_to = 0, copy_from = pop_size - 1;
	while (copy_to < copy_from) {
		while (num_replicas[copy_to] > 0)
			copy_to++;
		while (num_replicas[copy_from] <= 0)
			copy_from--;
		if (copy_to < copy_from) {
			pop_array[copy_to] = move(pop_array[copy_from]);
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
		for (j = 0; j < num_replicas[copy_from] - 1; j++) {
			pop_array[copy_to] = move(pop_array[copy_from]);
			copy_to++;
		}
		copy_from++;
	}
}

void Population::run(void)
{
    gsl_rng *r_thread[NUM_THREADS];
    FILE *fp = fopen("run.dat", "w");
    
	for (int i = 0; i < NUM_THREADS; i++) 
    {
        
		int tmp_rnd = gsl_rng_uniform_int(r, 10000000) + 1000000;
		initialize_rng(&r_thread[i], tmp_rnd);
	}
    
	double Beta = 0.0;
    double T = T_init;
	int num_sweeps;

    
    
    while (T != T_final) // Where the actual annealing happens
    { 
		fprintf(fp, "%5.5f\n", Beta);
		fflush(fp);
        num_sweeps = Beta <= 0.5 ? (SWEEPS / 4) : SWEEPS; // Do less sweeps for higher temperatures

        reSample(T, r);

        #pragma omp parallel for shared(pop_array, neighbor_table, Beta, r_thread, num_sweeps)
			for (int m = 0; m < pop_size; m++) 
            {
				int thread = omp_get_thread_num(); // BIG ISSUE HERE I THINK ... i think its ok now
				pop_array[m].doWolffAlgo(/*r_thread[thread], */neighbor_table, T, num_sweeps);
				pop_array[m].getTotalEnergy();
                cout << "Lattice " << m << " done running: Thread no. " << thread << "!";
            }   
    }   cout << "Done for T = " << T << "!";
}