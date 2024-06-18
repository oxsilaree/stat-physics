#include "population_class.h"

// Initializers
Population::Population(void) 
{
	throw invalid_argument("Invalid PA simulation initialization parameters.");
}

Population::Population(int nom_pop, int neighbor_table[L][L][nn_max][dim])
{
    // Initialize population
    Population::nom_pop = nom_pop;
    Population::max_pop = int(sqrt(nom_pop) * 10);
    Population::pop_size = nom_pop;
    Population::pop_array = unique_ptr<Lattice[]>(new Lattice[Population::max_pop]); // Not sure how this works but this makes
                                                                                     // the population
}