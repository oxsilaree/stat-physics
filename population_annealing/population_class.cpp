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

Population::Population(int nom_pop, gsl_rng *r, int nn_table[LEN][LEN][NN_MAX][DIM], double kappa)
{
    // Initialize population
    Population::kappa = kappa;
    Population::nom_pop = nom_pop;
    Population::max_pop = nom_pop + int(sqrt(nom_pop) * 10);
    Population::pop_size = nom_pop;
    Population::pop_array.reserve((int)sqrt(nom_pop)*10); // Reserve initial memory for the vector
    
    Population::energy_data.reserve((int)T_ITER);
    Population::magnetization_data.reserve((int)T_ITER);
    Population::spec_heat_data.reserve((int)T_ITER);
    Population::susceptibility_data.reserve((int)T_ITER);
    Population::temperature_values.reserve((int)T_ITER);

    for (int i = 0; i < nom_pop; ++i) {
        pop_array.push_back(Lattice(kappa));
    }
    // Population::pop_array = std::make_unique<Lattice[]>(Population::max_pop);
    // Population::pop_array = unique_ptr<Lattice[]>(new Lattice[Population::max_pop]); // Not sure how this works but this makes
                                                                                     // the population
    for (int i = 0; i < LEN; ++i) {
        for (int j = 0; j < LEN; ++j) {
            for (int pos = 0; pos < NN_MAX; ++pos) {
                for (int d = 0; d < DIM; ++d) {
                    Population::neighbor_table[i][j][pos][d] = nn_table[i][j][pos][d];
                }
            }
        }
    }

    Population::r = r;	    
    double dummy_1, dummy_2;
	energy_calcs(&dummy_1, &dummy_2);                
}

void Population::reSample(double T, gsl_rng *r, double avg_e, double var_e)
{
    int j = 0,     floor = 0, ceiling = 0;
    int new_pop_size = 0;
    double energy_j, weight_j, tau_j, expected_copies_j; // used for reweighting
    double Beta = 1/T;
    double T_prime = T - (T_INIT-T_FINAL)/T_ITER;
    if (T_prime < (T_INIT-T_FINAL)/T_ITER)
    {
        T_prime = 0.001;
    }
    // cout << "T_prime = " << T_prime << ".\n";
    double Beta_prime = 1/(T_prime);
    double d_Beta = Beta_prime - Beta; // Get the difference in previous and new inverse temperature
    //double d_Beta = CULLING_FRAC *sqrt(2 *PI /var_e); // This is a slightly more optimized way to run Pop.Annealing
	long double config_weight[pop_size];
	long double Q = 0.0;
	int num_replicas[pop_size];

    for (j = 0; j < pop_size; j++) // Get Q(beta, beta') for each lattice in population
    {
        Lattice* lattice_j = &pop_array[j];
        lattice_j->updateTotalEnergy(neighbor_table);
        energy_j = lattice_j->getTotalEnergy();
        weight_j = -d_Beta * (energy_j - avg_e); // These get massive
        config_weight[j] = exp(weight_j); // Storing these values is good to get tau later on
        Q += config_weight[j];
        cout << "Intermediate Q value: "<< Q << ". We just had energy " << energy_j << " and weight " << config_weight[j] << ".\n";
    }
    Q /= (double)pop_size;
    cout << "Q = " << Q << ". \n";

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
        // cout << "Number of replicas for " << j << "'th lattice: " << num_replicas[j] << ".\n";
    }
    try {
    if (new_pop_size > max_pop)
		throw "Maximum population size exceeded."; // Must 'catch' errors like this later on
    else if (new_pop_size <= 0)
        throw "Population size went to 0.";
    }
    catch (const char* errorMessage) {
        std::cerr << "Error: " << errorMessage << std::endl;
        exit(1);// Exit the program with an error status
    }
    // Vector routine for making a new population
    std::vector<Lattice> new_lattices;
    new_lattices.reserve(max_pop); // CHECK THIS LATER

    for (int i = 0; i < pop_size; i++) // Over set of replicas
    {
        for (int j = 0; j < num_replicas[i]; j++)
        {
            new_lattices.push_back(pop_array[i]);
        }
    }
    pop_array = new_lattices;
    pop_size = new_pop_size;
    cout << "Population size = " << pop_size << ".\n";
}

void Population::run(string kappastr)
{
    gsl_rng *r_thread[NUM_THREADS];
    double avg_e = 0.0, var_e = 0.0;
    FILE *fp = fopen("run.dat", "w");
    
	for (int i = 0; i < NUM_THREADS; i++) 
    {
		int tmp_rnd = gsl_rng_uniform_int(r, 10000000) + 1000000;
		initializeRNG(&r_thread[i], tmp_rnd);
	}
    
	double Beta = 0.0;
    double T = T_INIT;
	int num_sweeps;

    // Initialization
    for (int j = 0; j < nom_pop; j++)
    {
        Lattice* lattice_j = &pop_array[j];
        lattice_j->initializeSites();
    }
    
    while (T > T_FINAL) // Where the actual annealing happens
    { 
        Beta = 1/T;
        num_sweeps = Beta <= 0.5 ? (SWEEPS / 4) : SWEEPS; // Do less sweeps for higher temperatures

        collectData(T);
        energy_calcs(&avg_e, &var_e);
        reSample(T, r, avg_e, var_e);
        
        // Parallel version
        #pragma omp parallel for shared(pop_array, neighbor_table, Beta, r_thread, num_sweeps)
        for (int m = 0; m < pop_size; m++) 
        {
            int thread = omp_get_thread_num(); // check things regarding threads
            pop_array[m].doWolffAlgo(neighbor_table, T, num_sweeps);
            // pop_array[m].getTotalEnergy();
            // cout << "Lattice " << m << " done running: Thread no. " << thread << "!\n";
        }   

        
       
        cout << "Done for T = " << T << "!\n=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n";
        T -= double(T_INIT/T_ITER); // "Cooling" the system
        T = floor((100.*T)+.5)/100;
    } 
    // Load up the data into readable files  
    loadData(kappastr);
}

void Population::energy_calcs(double *avg_e, double *var_e) {
	int min_index = 0;
	double temp_min = 0.0;
	double avg_e2 = 0.0;
	double temp_avg_e = 0.0;
	for (int m = 0; m < pop_size; m++) {
        Lattice* lattice_m = &pop_array[m];
        lattice_m->updateTotalEnergy(neighbor_table);
        double temp_e = lattice_m->getTotalEnergy();
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
	// gs = pop_array[min_index];
	// gs_e = gs.get_energy();
}

void Population::collectData(double T)
{
    int m;
    double Beta = 1/T;
    double E,M,E_sq,M_sq,M_abs, C,X; // Placeholders
    double energy = 0, energy_sq = 0, mag = 0, mag_sq = 0, mag_abs = 0, spec_heat = 0, susc = 0;
    for (m = 0; m < pop_size; m++)
    {
        Lattice* lattice_m = &pop_array[m];
        lattice_m->updateTotalEnergy(neighbor_table);
        lattice_m->updateTotalMag();
        E = lattice_m->getTotalEnergy();
        M = lattice_m->getTotalMag();
        // E_sq = E*E;
        // M_sq = M*M;
        M_abs = abs(M);

        energy += E;
        energy_sq += E*E;
        mag += M;
        mag_sq += M*M;
        mag_abs += M_abs;
    }
    spec_heat = (energy_sq/pop_size - (energy/pop_size)*(energy/pop_size)) * (Beta*Beta);
    susc      = (mag_sq/pop_size - (mag_abs/pop_size)*(mag_abs/pop_size)) * Beta;
    temperature_values.push_back(T);
    energy_data.push_back(energy/pop_size);
    magnetization_data.push_back(mag_abs/pop_size); // Mag abs should give something nicer
    spec_heat_data.push_back(spec_heat);
    susceptibility_data.push_back(susc);
    // spec_heat_data.push_back(((energy_sq/pop_size) - pow((energy/pop_size),2)) * (Beta*Beta));
    // susceptibility_data.push_back(((mag_sq/pop_size) - pow((mag_abs/pop_size),2)) * Beta);
}

void Population::loadData(string kappastr)
{
    
    vector<double> T = temperature_values;
    vector<double> E = energy_data;
    vector<double> M = magnetization_data;
    vector<double> C = spec_heat_data;
    vector<double> X = susceptibility_data;

    ofstream emcx_data;
    emcx_data.open("data/emcx_data_" + kappastr + "_kappa.csv");
    emcx_data << "T,";
    for (vector<double>::iterator it=T.begin(); it != T.end(); ++it)
        {
            emcx_data << *it << ",";
        }
    emcx_data << "\n";
    emcx_data << "E,";
    for (vector<double>::iterator it=E.begin(); it != E.end(); ++it)
        {
            emcx_data << *it << ",";
        }        
    emcx_data << "\n";
    emcx_data << "M,";
    for (vector<double>::iterator it=M.begin(); it != M.end(); ++it)
        {
            emcx_data << *it << ",";
        }     
    emcx_data << "\n";
    emcx_data << "C,";
    for (vector<double>::iterator it=C.begin(); it != C.end(); ++it)
        {
            emcx_data << *it << ",";
        }       
    emcx_data << "\n";
    emcx_data << "X,";
    for (vector<double>::iterator it=X.begin(); it != X.end(); ++it)
        {
            emcx_data << *it << ",";
        }    
    emcx_data.close();
}