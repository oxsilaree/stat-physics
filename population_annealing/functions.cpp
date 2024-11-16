#include "functions.h"


int* getNeighbor(int i, int j, int pos)
{ // 0 -> left, 1 -> right, 2 -> up, 3 -> down, 4 -> up(NNN), 5 -> down(NNN)
    static int coords[2] = {0,0};
    if (pos == 0) { // left
        coords[1] = j;
        if (i == 0){
            coords[0] = i + (LEN-1);
        } else {
            coords[0] = i - 1;  
        }
    }
    if (pos == 1) { // right
        coords[1] = j;
        if (i == LEN - 1){
            coords[0] = i - (LEN-1);
        } else {
            coords[0] = i + 1;  
        }
    }
    if (pos == 2) { // up
        coords[0] = i;
        if (j == 0){
            coords[1] = j + (LEN-1);
        } else {
            coords[1] = j - 1;  
        }
    }
    if (pos == 3) { // down
        coords[0] = i;
        if (j == LEN - 1){
            coords[1] = j - (LEN-1);
        } else {
            coords[1] = j + 1;  
        }
    }
    if (pos == 4) { // up NNN
        coords[0] = i;
        if (j <= 1){
            coords[1] = j + (LEN-2);
        } else {
            coords[1] = j - 2;  
        }
    }
    if (pos == 5) { // down NNN
        coords[0] = i;
        if (LEN-2 <= j){
            coords[1] = j - (LEN-2);
        } else {
            coords[1] = j + 2;  
        }
    }
    return coords;
}

int neighbor_table[LEN][LEN][NN_MAX][DIM];

void makeNeighborTable()
{
    int* p;
    for (int i = 0; i < LEN; i++) // Make neighbor table
    {
        for (int j = 0; j < LEN; j++)
        {
            for (int pos = 0; pos < NN_MAX; pos++)
            {
                p = getNeighbor(i, j, pos); 
                neighbor_table[i][j][pos][0] = p[0];
                neighbor_table[i][j][pos][1] = p[1];
                // cout << neighbor_table[i][j][pos][0] << ", " << neighbor_table[i][j][pos][1] << endl;  
            }
        }
    }
}

void initializeRNG(gsl_rng **r, int seed) 
{
	gsl_rng_env_setup();
	const gsl_rng_type *T;
	T = gsl_rng_mt19937;
	*r = gsl_rng_alloc(T);
	gsl_rng_set(*r, seed);

	for (int i = 0; i < 10000; i++)
		gsl_rng_uniform(*r);
}




// Timer function
TimePoint timeCheck(TimePoint referenceTime) {
    TimePoint currentTime = std::chrono::steady_clock::now();
    chrono::duration<double> duration = currentTime - referenceTime;
    cout << fixed << std::setprecision(4);
    cout << "~~~~~~~~~~~~~~~~~~~~   Last portion took " << duration.count() << " seconds" << std::endl;
    
    return currentTime;
}

// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
const string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%Y-%m-%d_%H-%M-%S", &tstruct);

    return buf;
}

