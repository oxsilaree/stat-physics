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

Population::Population(int nom_pop, gsl_rng *r, double kappa, string mode)
{
    // Initialize population
    Population::kappa = kappa;
    Population::mode = mode; // Population or simulated annealing
    Population::nom_pop = nom_pop;
    Population::max_pop = nom_pop + int(sqrt(nom_pop) * 10);
    Population::pop_size = nom_pop;
    Population::avg_cluster_size = 0;
    Population::avg_nowrap_cluster_size = 0;
    Population::nowrap_count = 0;
    Population::r = r;
    Population::smoothed_var_e = 0;
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
        pop_array.push_back(Lattice(kappa,i));
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
    double dummy_1, dummy_2;
	calculateEnergies(&dummy_1, &dummy_2);                
}

void Population::reSample(double *Beta, gsl_rng *r, double avg_e, double var_e)
{
    int j = 0,     floor = 0, ceiling = 0;
    int new_pop_size = 0;
    double energy_j, weight_j, tau_j, expected_copies_j; // used for reweighting
    int new_family_counter = 0;

    double d_Beta = CULLING_FRAC * sqrt(2 *PI / var_e); // This is a slightly more optimized way to run Pop.Annealing
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
        weight_j = -d_Beta * (energy_j); // These get massive (-avg_e)
        config_weight[j] = exp(weight_j); // Storing these values is good to get tau later on
        Q += config_weight[j];
    }
    Q /= (double)pop_size;
    free_energy += log(Q);

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
    

  /*    *** Vector routine for changing population: do some reassignments and erases (this saves memory)
    Say we have population [A,B,C,D,E] with deleters = [B,D] and num_replicas = [2,0,3,0,1]. This code block will cause 
    the following progression: [A,B,C,D,E] --> [A,B,C,A,E] --> [A,C,C,A,E] --> [A,C,C,A,E,C] (population grows)

    Another demonstrative case is population [A,B,C,D,E] with deleters = [A,D] and num_replicas = [0,1,2,0,1].
    This will result in: [A,B,C,D,E] --> [A,B,C,C,E] --> [B,C,C,E] (population shrinks)

    The replacements and deletions occur from right to left due to deleters being a stack. The addition of
    a lattice is left->right because of push_back. This can lead to a relatively  mixed array of replicas.

    *//*
    for (int i = 0; i < pop_size; i++) // Iterate over current set of replicas
    {      
        pop_array[i].setNewFamily(new_family_counter);
        pop_array[i].updateRecentFamilies(new_family_counter);
        for (int j = 0; j < num_replicas[i] - 1; j++)   // Iterate over a single replica's number of new copies
        {                                               // minus 1 to ensure we don't redundantly add in the original replica
            if (deleters.empty() == false) 
            {                                           // deleters is a stack of indices for replicas to be deleted
                pop_array[deleters.top()] = pop_array[i];
                pop_array[deleters.top()].setNewFamily(new_family_counter); // Rewrite new family to ensure no identical replicas are paired later
                pop_array[deleters.top()].updateRecentFamilies(new_family_counter);
                deleters.pop();                         
            } else {
                pop_array.push_back(pop_array[i]);
                Lattice back_lattice = pop_array.back();
                back_lattice.setNewFamily(new_family_counter);      // BE CAREFUL THIS IS NEW CODE !!!!!!!!!!!!!!!!!
                back_lattice.updateRecentFamilies(new_family_counter);
                new_pop_size += 1;
            }
        }
        new_family_counter += 1;
    }
    while (deleters.empty() == false)
    {
        pop_array.erase(pop_array.begin() + deleters.top());
        deleters.pop();
        new_pop_size -= 1;
    }
    */

    // THIS IS A DIFFERENT LIFE/DEATH ROUTINE WHERE NEWBORNS ARE PLACED TO THE RIGHT OF THEIR ANCESTORS
    vector<Lattice>::iterator pop_it = pop_array.begin();  
    // int ii = 0;
    for (int ii = 0; ii < pop_size; ii++) // Iterate over current set of replicas
    {   
        pop_it->setNewFamily(new_family_counter);
        pop_it->updateRecentFamilies(new_family_counter);
        if (num_replicas[ii] == 0) {
            pop_array.erase(pop_it);
            new_pop_size -= 1;
        } else if (num_replicas[ii] == 1) {
            pop_it += 1;
        } else {
            for (int j = 0; j < num_replicas[ii] - 1; j++) {
                vector<Lattice>::iterator pop_next = pop_it + 1;
                pop_array.insert(pop_next, *pop_it);
                pop_next->setNewFamily(new_family_counter);
                pop_next->updateRecentFamilies(new_family_counter);
                new_pop_size += 1;
            }
        }
        new_family_counter += 1;
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

    for (int i = 0; i < NUM_THREADS; i++) {
		int tmp_rnd = gsl_rng_uniform_int(r, 10000000) + 1000000;
		initializeRNG(&r_thread[i], tmp_rnd);
	}
	
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
        int step_const;
        int num_sweeps;
        if (Beta == 0) {
            step_const = LEN*LEN;
        } else {
            step_const = ceil((LEN*LEN)/(sqrt(clustersize_data.back())));
        }

        num_sweeps = (Beta < 0.30 || Beta > 1.1) ? 1 : 20;
        num_steps = num_sweeps * step_const;

        double *in, *out;
        in = (double*) fftw_alloc_real(LEN);
        out = (double*) fftw_alloc_real(LEN);
        const fftw_plan p = fftw_plan_r2r_1d(LEN, in, out, FFTW_R2HC, FFTW_MEASURE);
        // omp_set_num_threads(NUM_THREADS);
        #pragma omp parallel for shared(pop_array, Beta, num_steps, p) schedule(auto)// , r_thread)
        for (int m = 0; m < pop_size; m++) 
        {
            pop_array[m].doWolffAlgo(&Beta, p, num_steps);
        }

        calculateEnergies(&avg_e, &var_e);
        calculateFamilies();
                
        // Cleanup
        fftw_destroy_plan(p);
        fftw_cleanup();
        
        measureOverlap();
        collectData(&Beta, avg_e, var_e);
        // makeHistograms(kappastr);
        reSample(&Beta, r, avg_e, var_e);
        /*
        if (var_e != 0)
        {
            collectData(&Beta, avg_e, var_e);
            reSample(&Beta, r, avg_e, var_e);
        } else {
            cout << "Variance went to 0. Sufficient data gathered. Beta = " << Beta << ".";
            cout << "\n=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n";
            break;
        }
        */
        
        
        
        
        cout << "Done for beta = " << Beta << "!\n=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n";

        // T -= double((T_INIT - T_FINAL)/T_ITER); // "Cooling" the system
        // T = floor((100.*T)+.5)/100;
    } 
    // Load up the data into readable files  
    loadData(kappastr);

    Lattice* lattice = &pop_array[0];
    lattice->printLattice();
}

void Population::calculateEnergies(double *avg_e, double *var_e) {
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
	Lattice* gs = &pop_array[min_index];
	gs_e = gs->getTotalEnergy();
}

void Population::collectData(double *Beta, double avg_e, double var_e)
{
    int m;
    double half_pop = pop_size/2;
    // double Beta = 1/T;
    double E,M, E2, M2, M_abs,C,X; // Placeholders
    double ene = 0, ene_sq = 0, mag = 0, mag_sq = 0, mag_abs = 0, spec_heat = 0, susc = 0;
    double CS, avg_clust_size = 0; // For wrapping
    double NWCS, avg_nowrap_clust_size = 0;
    double NWC, nowrap_counter = 0;
    double FR, AMP, freqs = 0, amps = 0;
    double OL, OA, OV;
    double ove, ove_abs, ove_var;
    for (m = 0; m < pop_size; m++)
    {
        Lattice* lattice_m = &pop_array[m];
        lattice_m->updateTotalEnergy();
        lattice_m->updateTotalMag();
        E = lattice_m->getTotalEnergy();
        M = lattice_m->getTotalMag();
        FR = lattice_m->getDomFreq();
        AMP = lattice_m->getDomAmplitude();
    
        M_abs = abs(M);
        ene = ene + E;
        ene_sq = ene_sq + (E*E);
        mag += M;
        mag_sq += M*M;
        mag_abs += M_abs;
        freqs += FR;
        amps += AMP;

        if (mode != "t") {
            CS = lattice_m->getAvgClusterSize();
            NWCS = lattice_m->getAvgNowrapClusterSize();
            NWC = lattice_m->getNoWrapCount();
            avg_clust_size += CS;
            avg_nowrap_clust_size += NWCS;
            nowrap_counter += NWC;
        }
    }

    if (mode == "t") { // I'm so sorry that these are named this way. LHS are for data collection only
            avg_clust_size = 2.0*avg_cluster_size;                     // RHS are data members of the population
            avg_nowrap_clust_size = 2.0*avg_nowrap_cluster_size;      // 2x to account for pop_size instead of half_pop in denom.
            nowrap_counter = nowrap_count;                         // when pushing to data file
    }
    double wrap_percent;
    if (mode == ("t"))
        wrap_percent = (wrap_counter)/(half_pop*num_steps);
        // wrap_percent = (double)((double)pop_size*num_steps/2 - (double)nowrap_counter)/((double)pop_size*num_steps/2);
    else
        wrap_percent = (double)(wrap_counter/(double)(pop_size*num_steps));

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
    clustersize_data.push_back(avg_clust_size/(pop_size*num_steps));
    nowrapclustersize_data.push_back(avg_nowrap_clust_size/(pop_size*num_steps));
    wrapping_data.push_back(wrap_percent);
    fft_freq_data.push_back(freqs/pop_size);
    fft_amp_data.push_back(amps/pop_size);
    rho_t_data.push_back(rho_t);
    unique_families_data.push_back(unique_families);
    overlap_data.push_back(overlap);
    abs_overlap_data.push_back(abs_overlap);
    var_overlap_data.push_back(var_overlap);
    free_energy_data.push_back(free_energy);

    cout << "Average (non-wrapping) cluster size: " << clustersize_data.back() << " (" << nowrapclustersize_data.back() << ").\n";
    cout << "E = " << ene/(pop_size*LEN*LEN) << ", C = " << spec_heat << ".\n";
    cout << "M = " << mag_abs/(pop_size*LEN*LEN) << ", X = " << susc << ".\n";
    cout << "Dom. Freq. = " << fft_freq_data.back() << ", Dom. Amplitude. = " << fft_amp_data.back() << "\n";
    cout << "Rho_T = " << rho_t << ", No. of Families = " << unique_families << ".\n";
    cout << "Absolute Overlap = " << abs_overlap << "\n";
    if (mode == "t") {
        cout << "Total number of steps (% w/ wrapping): " << pop_size*num_steps/2 << " (" << 100*wrap_percent <<"%).\n";
    } else {
        cout << "Total number of steps (% w/ wrapping): " << pop_size*num_steps << " (" << 100*wrap_percent <<"%).\n";
    }
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
    vector<double> RT = rho_t_data;
    vector<double> UF = unique_families_data;
    vector<double> OL = overlap_data;
    vector<double> OA = abs_overlap_data;
    vector<double> OV = var_overlap_data;
    vector<double> FE = free_energy_data;



    ofstream run_info;
    ofstream emcx_data;
    string lenstr = to_string(LEN);
    run_info.open("./data/" + mode + "/parameter_info_" + kappastr + "_kappa_" + lenstr + "_L" + ".csv"); // in SLURM
    // run_info.open("/Users/shanekeiser/Documents/ANNNI/populationannealing/data/" + mode + "/parameter_info_" + kappastr + "_kappa_" + lenstr + "_L" + ".csv"); // in my computer
    run_info << "Kappa,L,Initial Pop. Size,Culling Fraction\n";
    run_info << kappa << "," << LEN << "," << INIT_POP_SIZE << "," << CULLING_FRAC;
    run_info.close();

    emcx_data.open("./data/" + mode + "/emcx_data_" + kappastr + "_kappa_" + lenstr + "_L" + ".csv");
    // emcx_data.open("/Users/shanekeiser/Documents/ANNNI/populationannealing/data/" + mode + "/emcx_data_" + kappastr + "_kappa_" + lenstr + "_L" + ".csv");
    emcx_data << "Beta,Energy,Energy Squared,Magnetization,Magnetization Squared,Absolute Magnetization,";
    emcx_data << "Specific Heat,Susceptibility,Cluster Size,Non-Wrapping Cluster Size,Wrapping Probability,";
    emcx_data << "Dominant Frequency,Dominant Amplitude,Rho T,Unique Families,";
    emcx_data << "Overlap,Absolute Overlap,Overlap Variance,Free Energy\n";
    
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
    vector<double>::iterator it13 = RT.begin();
    vector<double>::iterator it14 = UF.begin();
    vector<double>::iterator it15 = OL.begin();
    vector<double>::iterator it16 = OA.begin();
    vector<double>::iterator it17 = OV.begin();
    vector<double>::iterator it18 = FE.begin();
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
                << *it12 << "," \
                << *it13 << "," \
                << *it14 << "," \
                << *it15 << "," \
                << *it16 << "," \
                << *it17 << "," \
                << *it18 << "\n";
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
        ++it13;
        ++it14;
        ++it15;
        ++it16;
        ++it17;
        ++it18;
    }     
    emcx_data.close();  
}

void Population::calculateFamilies(void) {
	int i;
	int family_size_dist[nom_pop];
	int family_size[nom_pop];

	for (i = 0; i < nom_pop; i++) {
		family_size[i] = 0;
		family_size_dist[i] = 0;
	}

//	Count the number of ancestors for each family
	for (i = 0; i < pop_size; i++) {
        Lattice* lattice_i = &pop_array[i];
		family_size[lattice_i->getFamily()]++;
    }
//	Count the distribution of family sizes
	for (i = 0; i < pop_size; i++) { // CHECK IF POP_SIZE OR NOM_POP: NOM_POP gives error if pop_size < nom_pop (naturally)
        Lattice* lattice_i = &pop_array[i];
		family_size_dist[lattice_i->getFamily()]++;
    }
	int counter = 0;
	rho_t = 0;

	for (i = 0; i < nom_pop; i++) {
		if (family_size[i] > 0) {
			counter++;
			rho_t += (double) family_size[i] *family_size[i] /(pop_size * pop_size);
		}
	}

	unique_families = counter;
	rho_t *= nom_pop;
}

int Population::countFamilies(void) {
	int i;
	int family_size[nom_pop];

	for (i = 0; i < nom_pop; i++) {
		family_size[i] = 0;
	}

//	Count the number of ancestors for each family
	for (i = 0; i < pop_size; i++) {
        Lattice* lattice_i = &pop_array[i];
		if (lattice_i->getTotalEnergy() == gs_e)
			family_size[lattice_i->getFamily()]++;
    }
	int counter = 0;
	for (i = 0; i < nom_pop; i++) {
		if (family_size[i] > 0)
			counter++;
	}

	return counter;
}


void Population::measureOverlap(){
    int r;
    double q, abs_q = 0;
    double overlap_sq = 0;
    overlap = 0, abs_overlap = 0, var_overlap = 0;
    for (int m = 0; m < pop_size; m++) {
        r = m;
        q = 0, abs_q = 0;
        Lattice* lattice_1 = &pop_array[m];

        while (r == m)
            r = rand() % pop_size;
        Lattice* lattice_2 = &pop_array[r];
    
        for (int i = 0; i < LEN; i++) {
            for (int j = 0; j < LEN; j++) {
                spinSite* spin_1 = lattice_1->getSpinSite(i,j);
                spinSite* spin_2 = lattice_2->getSpinSite(i,j);
                q += (spin_1->getSpin())*(spin_2->getSpin());
            }
        }
        q /= (LEN*LEN);
        abs_q = abs(q);
        overlap += q; // Definition of Overlap (C. Amey, J. Machta)
        abs_overlap += abs_q;
        overlap_sq += q*q;
    }
    overlap /= pop_size; // this is <q>
    abs_overlap /= pop_size; // this is <|q|>
    overlap_sq /= pop_size; // this is <q^2>
    var_overlap = overlap_sq - pow(abs_overlap,2); // this is <q^2> - <|q|>^2
}



void Population::runSA(string kappastr)
{
    double avg_e = 0.0, var_e = 0.0;

    double Beta = 0.0; // INITIAL BETA --> 'Infinite temperature' or some finite value


    // Data members

    // Initialization

    Lattice* lattice = &pop_array[0];
    lattice->initializeSites(&Beta);

    double *in, *out;
    in = (double*) fftw_alloc_real(LEN);
    out = (double*) fftw_alloc_real(LEN);
    const fftw_plan p = fftw_plan_r2r_1d(LEN, in, out, FFTW_R2HC, FFTW_MEASURE);
    
    while (Beta < MAX_BETA) // Where the actual annealing happens (T > T_FINAL), BETA_MAX = 5
    { 
        
        double E,M, E2, M2, M_abs,C,X; // Placeholders
        double ene = 0, ene_sq = 0, mag = 0, mag_sq = 0, mag_abs = 0, spec_heat = 0, susc = 0;
        double CS, avg_clust_size = 0; // For wrapping
        double NWCS, avg_nowrap_clust_size = 0;
        int NWC, nowrap_counter = 0;
        double FR, AMP, freqs = 0, amps = 0;
        
        int step_const;
        int num_sweeps;
        if (Beta == 0) {
            step_const = LEN*LEN;
        } else {
            step_const = ceil((LEN*LEN)/(clustersize_data.back()));
        }
        num_sweeps = (Beta < 0.30 || Beta > 1.1) ? 1 : 20;
        num_steps = num_sweeps * step_const;

        int blocks = INIT_POP_SIZE;

        for (int m = 0; m < blocks; m++) // Pop size is number of blocks here
        {
            pop_array[0].doWolffAlgo(&Beta, p, num_steps); /////// This is the actual algorithm being done
            
            Lattice* lattice = &pop_array[0];
            lattice->updateTotalEnergy();
            lattice->updateTotalMag();
            E = lattice->getTotalEnergy();
            M = lattice->getTotalMag();
            CS = lattice->getAvgClusterSize();
            NWCS = lattice->getAvgNowrapClusterSize();
            NWC = lattice->getNoWrapCount();
            FR = lattice->getDomFreq();
            AMP = lattice->getDomAmplitude();

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

        double wrap_percent = (double)(blocks*num_steps - nowrap_counter)/(blocks*num_steps);
    
        double var_e = (ene_sq/(blocks*LEN*LEN)) - pow(ene/(blocks*LEN),2);
        spec_heat = ((ene_sq/(blocks*LEN*LEN)) - pow(ene/(blocks*LEN),2)) * (Beta * Beta);
        susc      = ((mag_sq/(blocks*LEN*LEN)) - pow(mag_abs/(blocks*LEN),2)) * Beta;
        beta_values.push_back(Beta);
        // Quantities PER SPIN will have LEN*LEN in the denominator
        energy_data.push_back(ene/(blocks*LEN*LEN)); //
        energy_sq_data.push_back(ene_sq/(blocks*LEN*LEN));
        magnetization_data.push_back(mag/(blocks*LEN*LEN)); // Mag abs should give something nicer
        magnetization_sq_data.push_back(mag_sq/(blocks*LEN*LEN)); 
        magnetization_abs_data.push_back(mag_abs/(blocks*LEN*LEN));
        spec_heat_data.push_back(spec_heat);
        susceptibility_data.push_back(susc);
        clustersize_data.push_back(avg_clust_size/(blocks*num_steps));
        nowrapclustersize_data.push_back(avg_nowrap_clust_size/(blocks*num_steps));
        wrapping_data.push_back(wrap_percent);
        fft_freq_data.push_back(freqs/blocks);
        fft_amp_data.push_back(amps/blocks);
        rho_t_data.push_back(rho_t);
        unique_families_data.push_back(unique_families);

        cout << "Average (non-wrapping) cluster size: " << avg_clust_size/(blocks*num_steps) << " (" << avg_nowrap_clust_size/(blocks*num_steps) << ").\n";
        cout << "E = " << ene/(blocks*LEN*LEN) << ", C = " << spec_heat << ".\n";
        cout << "M = " << mag_abs/(blocks*LEN*LEN) << ", X = " << susc << ".\n";
        cout << "Dom. Freq. = " << freqs/blocks << ", Dom. Amplitude. = " << amps/blocks << "\n";
        // cout << "Rho_T = " << rho_t << ", No. of Families = " << unique_families << ".\n";
        cout << "Total number of steps (% w/ wrapping): " << blocks*num_steps << " (" << 100*wrap_percent <<"%).\n";
            
        cout << "Done for beta = " << Beta << "!\n=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n";

        double d_Beta = 0.005; // This is a slightly more optimized way to run Pop.Annealing
        if (Beta + d_Beta < MAX_BETA)
		    Beta += d_Beta;
	    else
		    Beta = MAX_BETA;
    } 
    // Cleanup

    fftw_destroy_plan(p);
    fftw_cleanup();
    // Load up the data into readable files  
    loadData(kappastr);

    lattice->printLattice();

}

void Population::doTwoReplica(double padd1, double padd2, int num_steps, gsl_rng *r, int index1, int index2) // Wolff-style two-replica cluster algorithm
{ 
    for (int s = 0; s < num_steps; s++)
        doTwoRepStep(padd1, padd2, r, index1, index2);
}

void Population::doTwoRepStep(double padd1, double padd2, gsl_rng *r, int index1, int index2)
{
    double cluster_size;
    int sp, i, j,   lx, ly,     oldspin, newspin,   oldspin_2, newspin_2;
    int current_x, current_y,   nn_i, nn_j; // Lattice indices
    int root_x, root_y,     coord_x, coord_y,   pos_update_x, pos_update_y; // Coordinates for wrapping criteria
    int old_x, old_y,       new_x, new_y;       // More coordinates for wrapping criteria 
    vector<pair<int, int> > stacker;
    stacker.reserve(LEN*LEN);
    stack<int> cluster_x;
    stack<int> cluster_y;
    bool wrapping_crit = false;
    double randnum;

    Lattice* lattice_1 = &pop_array[index1];             // From front of randomized indices
    Lattice* lattice_2 = &pop_array[index2];  // From rear of randomized indices
    i = gsl_rng_uniform_int(r, LEN);
    j = gsl_rng_uniform_int(r, LEN);
    
    spinSite* root_spin = lattice_1->getSpinSite(i,j); // Only need to check wrapping in 1 replica since we add
    root_spin->AddToCluster();                         // the same spins from both replicas to each cluster
    root_spin->Triangulate(0,0);
    cluster_x.push(i);
    cluster_y.push(j);
    stacker.push_back(make_pair(i, j));
    cluster_size = 1;

    sp = 1;
    while (sp)
    {   /* Pull a site off the stack from the nearest neighbours*/
        --sp;
        current_x = stacker.back().first;
        current_y = stacker.back().second;
        stacker.pop_back();
        spinSite* current_spin = lattice_1->getSpinSite(current_x, current_y);
        spinSite* current_spin_2 = lattice_2->getSpinSite(current_x, current_y);
        oldspin = current_spin->getSpin(), oldspin_2 = current_spin_2->getSpin(); // Main difference is storing
        newspin = -oldspin, newspin_2 = -oldspin_2;                               // Spin values for both reps
        coord_x = current_spin->getX();
        coord_y = current_spin->getY();

        /* Check the neighbours using neighbour table */

        for (int k = 0; k < NN_MAX; k++)
        {   
            
            nn_i = neighbor_table[current_x][current_y][k][0];
            nn_j = neighbor_table[current_x][current_y][k][1];

            spinSite* neighbor_spin   = lattice_1->getSpinSite(nn_i,nn_j); 
            spinSite* neighbor_spin_2 = lattice_2->getSpinSite(nn_i,nn_j);
            bool neighbor_checked = neighbor_spin->checkStatus();
            if (wrapping_crit == false && neighbor_checked == true) { // Check for wrapping once

                old_x = neighbor_spin->getX();
                old_y = neighbor_spin->getY();
                int pos_update_x = coord_x, pos_update_y = coord_y;
                if (k == 0) pos_update_x -= 1;
                else if (k == 1) pos_update_x += 1;
                else if (k == 2) pos_update_y += 1;
                else if (k == 3) pos_update_y -= 1;
                else if (k == 4) pos_update_y += 2;
                else if (k == 5) pos_update_y -= 2;

                neighbor_spin->Triangulate(pos_update_x, pos_update_y);
                new_x = neighbor_spin->getX();
                new_y = neighbor_spin->getY();

                if (old_x != new_x || old_y != new_y) {
                    wrapping_crit = true;
                    wrap_counter++;
                }
            }   else if (neighbor_checked == false)  {
                //////// CHANGE CRITERION HERE (changed)
                randnum = gsl_rng_uniform(r);
                if ((neighbor_spin->getSpin() == oldspin && neighbor_spin_2->getSpin() == oldspin_2 && k <= 3 && randnum <= padd1) ||
                    (neighbor_spin->getSpin() == newspin && neighbor_spin_2->getSpin() == newspin_2 && k >= 4 && randnum <= padd2)) 
                {
                    int pos_update_x = coord_x, pos_update_y = coord_y;
                    if (k == 0) pos_update_x -= 1;
                    else if (k == 1) pos_update_x += 1;
                    else if (k == 2) pos_update_y += 1;
                    else if (k == 3) pos_update_y -= 1;
                    else if (k == 4) pos_update_y += 2;
                    else if (k == 5) pos_update_y -= 2;

                    sp += 1;
                    stacker.push_back(make_pair(nn_i, nn_j));
                    neighbor_spin->AddToCluster();
                    neighbor_spin->Triangulate(pos_update_x, pos_update_y);
                    cluster_x.push(nn_i);
                    cluster_y.push(nn_j);
                    cluster_size += 1;
                }
            }
        }
    }

    // For simpler percolation data, we are just interested in the non-wrapping cluster size.
    if (wrapping_crit == false)
    {
        avg_nowrap_cluster_size += cluster_size;
        nowrap_count++;
    }
    avg_cluster_size += cluster_size;
    /* Go over all the spins in the stack and flip if they are in the cluster */
    while (!cluster_x.empty()) {
        lx = cluster_x.top();
        ly = cluster_y.top();
        spinSite* spin_1 = lattice_1->getSpinSite(lx,ly);
        spinSite* spin_2 = lattice_2->getSpinSite(lx,ly);
        spin_1->Flip();
        spin_1->Reset();
        spin_2->Flip();
        spin_2->Reset();
        cluster_x.pop();
        cluster_y.pop();
    }
}

void Population::runTR(string kappastr)
{   
    gsl_rng *r_thread[NUM_THREADS];
    double avg_e = 0.0, var_e = 0.0;
    int tmp_rnd;
    int temp_steps = 0;
    auto ref_time = chrono::high_resolution_clock::now(); // PROFILER SET-UP
    
	for (int i = 0; i < NUM_THREADS; i++) 
    {
		tmp_rnd = gsl_rng_uniform_int(r, 10000000) + 1000000;
		initializeRNG(&r_thread[i], tmp_rnd);
	}
    
    double Beta = 0.0; // INITIAL BETA --> 'Infinite temperature' or some finite value
    // Initialization
    for (int j = 0; j < nom_pop; j++)
    {
        Lattice* lattice_j = &pop_array[j];
        lattice_j->initializeSites(&Beta);
    }
    
    while (Beta < MAX_BETA) // Where the actual annealing happens (T > T_FINAL), BETA_MAX = 5
    { 
        int step_const, num_sweeps;
        temp_steps++;
        if (clustersize_data.empty()) {
            step_const = LEN*LEN;
        } else if  (ceil((LEN*LEN)/(clustersize_data.back())) > 10) {
            step_const = ceil((LEN*LEN)/(clustersize_data.back()));
        } else {
            step_const = 10;
        }
        num_sweeps = (Beta < 0.35) ? 1 : (Beta < 1.1) ? 20 : 5;
        num_steps = num_sweeps * step_const;
        double *in, *out;
        in = (double*) fftw_alloc_real(LEN);
        out = (double*) fftw_alloc_real(LEN);
        const fftw_plan p = fftw_plan_r2r_1d(LEN, in, out, FFTW_R2HC, FFTW_MEASURE);
        // omp_set_num_threads(NUM_THREADS);
        
        // THIS IS THE ACTUAL RUN
        wrap_counter = 0, nowrap_count = 0;
        avg_cluster_size = 0, avg_nowrap_cluster_size = 0;

        // If odd number of replicas, remove a random replica
        if (pop_size % 2 == 1) {
            gsl_rng *r1 = gsl_rng_alloc(gsl_rng_mt19937); // You can replace gsl_rng_mt19937 with another RNG algorithm if desired
            gsl_rng_set(r1, time(NULL));
            int value = gsl_rng_uniform_int(r1, pop_size);
            pop_array.erase(pop_array.begin() + value);
            gsl_rng_free(r1);
            pop_size -= 1;
        }
        ref_time = chrono::high_resolution_clock::now(); // PROFILER STEP
        cout << "~~~~~~~~~~~~~~~~~~~~ Pairing and swapping...\n"; // PROFILER STEP

        // Pre-make random pairs of indices
        vector<int> indices(pop_size);
        for (int i = 0; i < pop_size; i++) {
            indices[i] = i;
        }   
        random_device rd;
        mt19937 g(rd());
        shuffle(indices.begin(), indices.end(), g);
        int half_pop = pop_size/2;

        
        // Swap replicas if we have correlated pairs
        
        if (Beta > 0.1) {
            for (int m = 0; m < half_pop; m++) {
                vector<int>::iterator it, it2, it3; // Use this to edit the indices vector
                int it_count = 2;
                it = indices.begin() + 2*m;                         // v--(d = 50)--v
                
                while (abs(indices[2*m] - indices[2*m + 1]) < MIN_DISTANCE) { // "If randomly paired replicas are within d replicas 
                    vector<int>::iterator swapper = indices.begin() + ((2*m + it_count) % pop_size); // from each other on the population array"
                    iter_swap(it+1, swapper);                       
                    it_count++;
                } 
                /*
                // This one is for k-step look-back
                while (haveSharedFamily(&pop_array[indices[2*m]], &pop_array[indices[2*m + 1]])) { // "If randomly paired replicas share a family
                    iter_swap(it+1, it+it_count);                                                       // up to k annealing steps back"
                    it_count++;
                }
                *//*
                // This one is for single look-back
                while (pop_array[indices[2*m]].getNewFamily() == pop_array[indices[2*m + 1]].getNewFamily()) { // "If randomly paired replicas are identical"
                    iter_swap(it+1, it+it_count);
                    it_count++;
                }
                */
                
            } 
        }
        
        
        ref_time = timeCheck(ref_time); // PROFILER STEP
        
        cout << "~~~~~~~~~~~~~~~~~~~~ Doing Wolff steps... \n";
        // Decorrelate pairs using Wolff steps
        
        double w_padd1 = 1 - exp(-2*Beta*J);
        double w_padd2 = 1 - exp(-2*Beta*J*kappa);
        int wolff_steps;
        if (temp_steps > 2) {
            wolff_steps = (spec_heat_data.back()*10) + 1;
            wolff_steps = max(30,wolff_steps);
            cout << "Number of wolff Steps = " << wolff_steps << "\n";
        } else {
            wolff_steps = 1;
        }
        #pragma omp parallel for shared(pop_array, w_padd1, w_padd2) schedule(dynamic, 20)
        for (int m = 0; m < pop_size; m++) {
            Lattice* lattice_m = &pop_array[m];
            for (int j = 0; j < 30; j++) {
                lattice_m->doStep(&w_padd1, &w_padd2);
            }
        }
        
        ref_time = timeCheck(ref_time); // PROFILER STEP
        
        double padd1 = 1 - exp(-4 * Beta * J);
        double padd2 = 1 - exp(-4 * Beta * J * kappa);
        
        cout << "~~~~~~~~~~~~~~~~~~~~ Doing Two-Replica steps...\n";
        // Pair up lattices and do two replica cluster moves (parallelized) --- only pair once
        
        #pragma omp parallel for shared(pop_array, padd1, padd2, num_steps, r_thread, indices, wrap_counter) schedule(dynamic, 10)
        for (int m = 0; m < half_pop; m++) {
            int thread = omp_get_thread_num();
            doTwoReplica(padd1, padd2, num_steps, r_thread[thread % NUM_THREADS], indices[2*m], indices[2*m + 1]);
        }                                            // For some reason, need to mod this above.
        ref_time = timeCheck(ref_time); // PROFILER STEP
        
        
        // Do two replica cluster moves but re-pair replicas every sweep.
        /*
        for (int sw = 0; sw < num_sweeps; sw++) {
            // Shuffle, swap and fix pairs once per sweep
            if (sw % 5 == 0)
                shuffle(indices.begin(), indices.end(), g);
            for (int m = 0; m < half_pop; m++) {
                vector<int>::iterator it; // Use this to edit the indices vector
                int it_count = 2;
                it = indices.begin() + 2*m;                         // v--(d = 50)--v
                while (abs(indices[2*m] - indices[2*m + 1]) < MIN_DISTANCE) { // "If randomly paired replicas are within d replicas 
                    vector<int>::iterator swapper = indices.begin() + (((2*m) + it_count) % pop_size); // from each other on the population array"
                    iter_swap(it+1, swapper);                       
                    it_count++;
                }
            }
            #pragma omp parallel for shared(pop_array, padd1, padd2, r_thread, indices, wrap_counter) schedule(dynamic, 10)
            for (int m = 0; m < half_pop; m++) {
                int thread = omp_get_thread_num();
                for (int st = 0; st < step_const; st++) {
                    doTwoRepStep(padd1, padd2, r_thread[thread % NUM_THREADS], indices[2*m], indices[2*m + 1]);
                }
            }
        }
        
        
        ref_time = timeCheck(ref_time); // PROFILER STEP
        */
        cout << "~~~~~~~~~~~~~~~~~~~~ Doing FFTs...\n";
        #pragma omp parallel for shared (p) schedule(dynamic, 10)
        for (int m = 0; m < pop_size; m++) {
            pop_array[m].doFFT(p);
        }
        
        // FFT plan Cleanup
        fftw_destroy_plan(p);
        fftw_cleanup();
        ref_time = timeCheck(ref_time); // PROFILER STEP

        cout << "~~~~~~~~~~~~~~~~~~~~ Doing data collection...\n";
        calculateEnergies(&avg_e, &var_e);
        if (smoothed_var_e == 0) {
            smoothed_var_e += var_e;
        } else {
            smoothed_var_e += var_e;
            smoothed_var_e /= 2;
        }
        // CHECK FOR ENERGY OUTLIERS 
        int anomaly_counter = 0;
        for (int mm = 0; mm < pop_size; mm++) {
            Lattice* lattice_mm = &pop_array[indices[mm]];
            double ene_mm = lattice_mm->getTotalEnergy();
            if (abs(ene_mm - avg_e) > 4*sqrt(var_e)){ // Find lattices that are 5 sigma away.
                // cout << "Lattice number " << indices[mm] << " is weird. " << mm << " is its position in the randomized order. E = " << ene_mm << ". Family number = " << lattice_mm->getFamily() << ", latest NewFamily number = " << lattice_mm->getNewFamily() << ".\n";
                anomaly_counter++;
            }
        }

        calculateFamilies();
        measureOverlap();
        if (Beta > 0) {
            collectData(&Beta, avg_e, var_e);
        }
        ref_time = timeCheck(ref_time); // PROFILER STEP
        // makeHistograms(kappastr);
        cout << "~~~~~~~~~~~~~~~~~~~~ Doing resampling...\n";
        reSample(&Beta, r, avg_e, var_e);
        ref_time = timeCheck(ref_time); // PROFILER STEP
        /*
        if (var_e == 0)
        {
            cout << "Variance went to 0. Sufficient data gathered. Beta = " << Beta << ".";
            cout << "\n=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n";
            break;
        }
        */
        
        
        
        
        cout << "Done for beta = " << Beta << "! StDev in energy = " << sqrt(var_e) << ", no. of anomalies (>4sigma) = " << anomaly_counter;
        cout << "\nNo. of steps = " << step_const << ", No. of sweeps = " << num_sweeps << "\n=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n";

        // T -= double((T_INIT - T_FINAL)/T_ITER); // "Cooling" the system
        // T = floor((100.*T)+.5)/100;
    } 
    // Load up the data into readable files  
    loadData(kappastr);

    Lattice* lattice = &pop_array[0];
    lattice->printLattice();
    cout << "Total temperature steps: " << temp_steps << "\n";
}

void Population::makeHistograms(string kappastr) {
    ofstream ene_histograms;
    ofstream mag_histograms;
    string lenstr = to_string(LEN);
    // ene_histograms.open("./data/" + mode + "/ene_hist" + kappastr + "_kappa_" + lenstr + "_L" + ".csv", ios::app); // in SLURM
    // ene_histograms.open("/Users/shanekeiser/Documents/ANNNI/populationannealing/data/" + mode + "/ene_hist_" + kappastr + "_kappa_" + lenstr + "_L" + ".csv",ios::app); // in my computer
    // mag_histograms.open("./data/" + mode + "/mag_hist_" + kappastr + "_kappa_" + lenstr + "_L" + ".csv",ios::app); // in SLURM
    // mag_histograms.open("/Users/shanekeiser/Documents/ANNNI/populationannealing/data/" + mode + "/mag_hist_" + kappastr + "_kappa_" + lenstr + "_L" + ".csv",ios::app); // in my computer

    for (int m = 0; m < max_pop; m++) {
        if (m < pop_size) {
            Lattice* lattice_m = &pop_array[m];
            ene_histograms << lattice_m->getTotalEnergy() << ",";
            mag_histograms << lattice_m->getTotalMag() << ",";
        } else {
            ene_histograms << "5000" << ",";
            mag_histograms << "5000" << ",";
        }
    }
    ene_histograms << "\n";
    mag_histograms << "\n";
    ene_histograms.close();
    mag_histograms.close();
}


// Function to check if two lattices share at least one family number
bool Population::haveSharedFamily(Lattice* lattice1, Lattice* lattice2) {
    // Create a set from lattice1's recent_families
    deque<int> rec_fams_1 = lattice1->getRecentFamilies();
    deque<int> rec_fams_2 = lattice2->getRecentFamilies();
    set<int> families_set(rec_fams_1.begin(), rec_fams_1.end());

    // Check if any element in lattice2's recent_families exists in the set
    for (int family : rec_fams_2) {
        if (families_set.count(family)) {
            return true;  // At least one shared number found
        } else {
            return false;
        }
    }
    return false;
}



