#include "parameters.h"
#include "spin_class.h"

void spinSite::AssignValues()
{
    spin = rand() % 2 == 0 ? 1 : -1; //: Assign a random value to spin (1 or -1)
    in_cluster = false; // Start off "not in cluster"
}

