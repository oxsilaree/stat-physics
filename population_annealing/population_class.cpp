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

Population::Population(int nom_pop, gsl_rng *r, double kappa)
{
    // Initialize population
    Population::kappa = kappa;
    Population::nom_pop = nom_pop;
    Population::max_pop = nom_pop + int(sqrt(nom_pop) * 10);
    Population::pop_size = nom_pop;
    Population::pop_array.reserve(max_pop); // Reserve initial memory for the vector
    
    Population::energy_data.reserve((int)2*T_ITER);
    Population::magnetization_data.reserve((int)2*T_ITER);
    Population::spec_heat_data.reserve((int)2*T_ITER);
    Population::susceptibility_data.reserve((int)2*T_ITER);
    Population::clustersize_data.reserve((int)2*T_ITER);
    Population::nowrapclustersize_data.reserve((int)2*T_ITER);
    Population::beta_values.reserve((int)2*T_ITER);
    Population::energy_sq_data.reserve((int)2*T_ITER);
    Population::magnetization_sq_data.reserve((int)2*T_ITER);
    Population::magnetization_abs_data.reserve((int)2*T_ITER);
    Population::wrapping_data.reserve((int)2*T_ITER);
    

    for (int i = 0; i < nom_pop; ++i) {
        pop_array.push_back(Lattice(kappa));
    }
    /*
    for (int i = 0; i < LEN; ++i) {
        for (int j = 0; j < LEN; ++j) {
            for (int pos = 0; pos < NN_MAX; ++pos) {
                for (int d = 0; d < DIM; ++d) {
                    Population::neighbor_table[i][j][pos][d] = nn_table[i][j][pos][d];
                }
            }
        }
    }
    */
    Population::r = r;	    
    double dummy_1, dummy_2;
	energy_calcs(&dummy_1, &dummy_2);                
}

void Population::reSample(double *Beta, gsl_rng *r, double avg_e, double var_e)
{
    int j = 0,     floor = 0, ceiling = 0;
    int new_pop_size = 0;
    double energy_j, weight_j, tau_j, expected_copies_j; // used for reweighting


    double d_Beta = CULLING_FRAC * LEN * sqrt(2 *PI /var_e); // This is a slightly more optimized way to run Pop.Annealing
    if (*Beta + d_Beta < MAX_BETA)
		*Beta += d_Beta;
	else
		*Beta = MAX_BETA;

	long double config_weight[pop_size];
	long double Q = 0.0;
	int num_replicas[pop_size];
    stack<int> deleters;

    for (j = 0; j < pop_size; j++) // Get Q(beta, beta') for each lattice in population
    {
        Lattice* lattice_j = &pop_array[j];
        lattice_j->updateTotalEnergy();
        energy_j = lattice_j->getTotalEnergy();
        weight_j = -d_Beta * (energy_j - avg_e); // These get massive
        config_weight[j] = exp(weight_j); // Storing these values is good to get tau later on
        Q += config_weight[j];
        // cout << "Intermediate Q value: "<< Q << ". We just had energy " << energy_j << " and weight " << weight_j << ".\n";
    }
    Q /= (double)pop_size;
    // cout << "Q = " << Q << ". \n";

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
            
        if (num_replicas[j] == 0)
        {
            deleters.push(j);
        }
        // new_pop_size += num_replicas[j]; // AND THISSSSSSSSS
        // cout << "Number of replicas for " << j << "'th lattice: " << num_replicas[j] << ".\n";
    }
    new_pop_size = pop_size; ////////// MAKE SURE TO CHANGE THIS BACK
    
/*
    // Vector routine for making a new population: just make a copy
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
*/
  
  // Vector routine for changing population: do some reassignments and erases (saves memory)
  /*
    Say we have population [A,B,C,D,E] with deleters = [B,D] and num_replicas = [2,0,3,0,1]. This code block will cause 
    the following progression: [A,B,C,D,E] --> [A,B,C,A,E] --> [A,C,C,A,E] --> [A,C,C,A,E,C] (population grows)

    Another demonstrative case is population [A,B,C,D,E] with deleters = [A,D] and num_replicas = [0,1,2,0,1].
    This will result in: [A,B,C,D,E] --> [A,B,C,C,E] --> [B,C,C,E] (population shrinks)

    The replacements and deletions occur from right to left due to deleters being a stack. 
    The addition of a lattice is from left to right because of the use of push_back.
    This can lead to a pretty mixed population.
  */
    for (int i = 0; i < pop_size; i++) // Iterate over current set of replicas
    {   
        for (int j = 0; j < num_replicas[i] - 1; j++)   // Iterate over a single replica's number of new copies
        {                                               // minus 1 to ensure we don't redundantly add in the original replica
            if (deleters.empty() == false) 
            {                                           // deleters is a stack of indices for replicas to be deleted
                pop_array[deleters.top()] = pop_array[i];
                deleters.pop();                         
            } else {
                pop_array.push_back(pop_array[i]);
                new_pop_size += 1;
            }
        }
    }
    while (deleters.empty() == false)
    {
        pop_array.erase(pop_array.begin() + deleters.top());
        deleters.pop();
        new_pop_size -= 1;
    }

    pop_size = new_pop_size;

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

    cout << "Population size = " << pop_size << ".\n";
}

void Population::run(string kappastr)
{
    gsl_rng *r_thread[NUM_THREADS];
    double avg_e = 0.0, var_e = 0.0;
    // FILE *fp = fopen("run.dat", "w");
    /*
	for (int i = 0; i < NUM_THREADS; i++) 
    {
		int tmp_rnd = gsl_rng_uniform_int(r, 10000000) + 1000000;
		initializeRNG(&r_thread[i], tmp_rnd);
	}
    */
	
    // double T = T_INIT;
    double Beta = 0.0; // INITIAL BETA --> 'Infinite temperature' or some finite value


    // Initialization
    for (int j = 0; j < nom_pop; j++)
    {
        Lattice* lattice_j = &pop_array[j];
        lattice_j->initializeSites(&Beta);
    }
    
    while (Beta < MAX_BETA) // Where the actual annealing happens (T > T_FINAL), BETA_MAX = 5
    { 
        // Parallel version
        // int num_sweeps = (Beta < .46) ? SWEEPS : SWEEPS/4;
        // int num_sweeps = (Beta < .25) ? SWEEPS*10 : (Beta < .46) ? SWEEPS : SWEEPS/4;
        int num_sweeps = 0;
        if (Beta < 0.35) {
            num_sweeps = SWEEPS;
        } else if (Beta < 0.7) {
            num_sweeps = SWEEPS * 3;
        } else {
            num_sweeps = SWEEPS / 4;
        }

        double *in, *out;
        in = (double*) fftw_alloc_real(LEN);
        out = (double*) fftw_alloc_real(LEN);
        const fftw_plan p = fftw_plan_r2r_1d(LEN, in, out, FFTW_R2HC, FFTW_MEASURE);
        #pragma omp parallel for shared(pop_array, Beta, num_sweeps, p)// , r_thread)
        for (int m = 0; m < pop_size; m++) 
        {
            // int thread = omp_get_thread_num(); // check things regarding threads
            pop_array[m].doWolffAlgo(&Beta, p);
            // pop_array[m].getTotalEnergy();
            // cout << "Lattice " << m << " done running: Thread no. " << thread << "!\n";
        }   
        
        energy_calcs(&avg_e, &var_e);
                
        // Cleanup
        fftw_destroy_plan(p);
        fftw_cleanup();
        
        if (var_e != 0)
        {
            collectData(&Beta, avg_e, var_e);
            reSample(&Beta, r, avg_e, var_e);
        } else {
            cout << "Variance went to 0. Sufficient data gathered. Beta = " << Beta << ".";
            cout << "\n=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n";
            break;
        }
        
        

        
        
        cout << "Done for beta = " << Beta << "!\n=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n";

        // T -= double((T_INIT - T_FINAL)/T_ITER); // "Cooling" the system
        // T = floor((100.*T)+.5)/100;
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
        lattice_m->updateTotalEnergy();
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

void Population::collectData(double *Beta, double avg_e, double var_e)
{
    int m;
    // double Beta = 1/T;
    double E,M, E2, M2, M_abs,C,X; // Placeholders
    double ene = 0, ene_sq = 0, mag = 0, mag_sq = 0, mag_abs = 0, spec_heat = 0, susc = 0;
    double CS, avg_clust_size = 0; // For wrapping
    double NWCS, avg_nowrap_clust_size = 0;
    int NWC, nowrap_counter = 0;
    double FR, AMP, freqs = 0, amps = 0;
    //int num_sweeps = (*Beta < 0.46) ? SWEEPS : SWEEPS/4;
    int num_sweeps = (*Beta < .25) ? SWEEPS*10 : (* Beta < .46) ? SWEEPS : SWEEPS/4;
    for (m = 0; m < pop_size; m++)
    {
        Lattice* lattice_m = &pop_array[m];
        lattice_m->updateTotalEnergy();
        lattice_m->updateTotalMag();
        E = lattice_m->getTotalEnergy();
        M = lattice_m->getTotalMag();
        CS = lattice_m->getAvgClusterSize();
        NWCS = lattice_m->getAvgNowrapClusterSize();
        NWC = lattice_m->getNoWrapCount();
        FR = lattice_m->getDomFreq();
        AMP = lattice_m->getDomAmplitude();
    
        // E_sq = E*E;
        // M_sq = M*M;
        M_abs = abs(M);
        ene = ene + E;
        ene_sq = ene_sq + (E*E);
        mag += M;
        mag_sq += M*M;
        mag_abs += M_abs;
        avg_clust_size += CS;
        avg_nowrap_clust_size += NWCS;
        nowrap_counter += NWC;
        freqs += FR;
        amps += AMP;
    }
    
    double wrap_percent = (double)(pop_size*num_sweeps*STEPS - nowrap_counter)/(pop_size*num_sweeps*STEPS);
    

    spec_heat = ((ene_sq/(pop_size*LEN*LEN)) - pow(ene/(pop_size*LEN),2)) * (*Beta * *Beta);
    susc      = ((mag_sq/(pop_size*LEN*LEN)) - pow(mag_abs/(pop_size*LEN),2)) * *Beta;
    beta_values.push_back(*Beta);
    // Quantities PER SPIN will have LEN*LEN in the denominator
    energy_data.push_back(ene/(pop_size*LEN*LEN)); //
    energy_sq_data.push_back(ene_sq/(pop_size*LEN*LEN));
    magnetization_data.push_back(mag/(pop_size*LEN*LEN)); // Mag abs should give something nicer
    magnetization_sq_data.push_back(mag_sq/(pop_size*LEN*LEN)); 
    magnetization_abs_data.push_back(mag_abs/(pop_size*LEN*LEN));
    spec_heat_data.push_back(spec_heat);
    susceptibility_data.push_back(susc);
    clustersize_data.push_back(avg_clust_size/(pop_size*num_sweeps*STEPS));
    nowrapclustersize_data.push_back(avg_nowrap_clust_size/(pop_size*num_sweeps*STEPS));
    wrapping_data.push_back(wrap_percent);
    fft_freq_data.push_back(freqs/pop_size);
    fft_amp_data.push_back(amps/pop_size);

    cout << "Average (non-wrapping) cluster size: " << avg_clust_size/(pop_size*num_sweeps*STEPS) << " (" << avg_nowrap_clust_size/(pop_size*num_sweeps*STEPS) << ").\n";
    cout << "E = " << ene/(pop_size*LEN*LEN) << ", C = " << spec_heat << ".\n";
    cout << "M = " << mag_abs/(pop_size*LEN*LEN) << ", X = " << susc << ".\n";
    cout << "Dom. Freq. = " << freqs/pop_size << ", Dom. Amplitude. = " << amps/pop_size << "\n";
    cout << "Total number of steps (% w/ wrapping): " << pop_size*num_sweeps*STEPS << " (" << 100*wrap_percent <<"%).\n";
}

void Population::loadData(string kappastr)
{
    vector<double> B = beta_values;
    vector<double> E = energy_data;
    vector<double> M = magnetization_data;
    vector<double> E2 = energy_sq_data;
    vector<double> MA = magnetization_abs_data;
    vector<double> M2 = magnetization_sq_data;
    vector<double> C = spec_heat_data;
    vector<double> X = susceptibility_data;
    vector<double> CS = clustersize_data;
    vector<double> NWCS = nowrapclustersize_data;
    vector<double> W = wrapping_data;
    vector<double> FR = fft_freq_data;
    vector<double> AM = fft_amp_data;


    ofstream run_info;
    ofstream emcx_data;
    string lenstr = to_string(LEN);
    run_info.open("data/parameter_info_" + kappastr + "_kappa_" + lenstr + "_L" + ".csv");
    run_info << kappa << "," << LEN << "," << INIT_POP_SIZE << "," << CULLING_FRAC << ",";
    run_info.close();

    emcx_data.open("data/emcx_data_" + kappastr + "_kappa_" + lenstr + "_L" + ".csv");
    emcx_data << "Beta,Energy,Energy Squared,Magnetization,Magnetization Squared,Absolute Magnetization,";
    emcx_data << "Specific Heat,Susceptibility,Cluster Size,Non-Wrapping Cluster Size,Wrapping Probability,";
    emcx_data << "Dominant Frequency, Dominant Amplitude\n";
    
    vector<double>::iterator it1 = E.begin();
    vector<double>::iterator it2 = E2.begin();
    vector<double>::iterator it3 = M.begin();
    vector<double>::iterator it4 = M2.begin();
    vector<double>::iterator it5 = MA.begin();
    vector<double>::iterator it6 = C.begin();
    vector<double>::iterator it7 = X.begin();
    vector<double>::iterator it8 = CS.begin();
    vector<double>::iterator it9 = NWCS.begin();
    vector<double>::iterator it10 = W.begin();
    vector<double>::iterator it11 = FR.begin();
    vector<double>::iterator it12 = AM.begin();
    for (vector<double>::iterator it=B.begin(); it != B.end(); ++it)      
    {
        emcx_data << *it << "," \
                << *it1 << "," \
                << *it2 << "," \
                << *it3 << "," \
                << *it4 << "," \
                << *it5 << "," \
                << *it6 << "," \
                << *it7 << "," \
                << *it8 << "," \
                << *it9 << "," \
                << *it10 << "," \
                << *it11 << "," \
                << *it12 << "\n";
        ++it1;
        ++it2;
        ++it3;
        ++it4;
        ++it5;
        ++it6;
        ++it7;
        ++it8;
        ++it9;
        ++it10;
        ++it11;
        ++it12;
    }     
    emcx_data.close();  
}