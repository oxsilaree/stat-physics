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

void initialize_rng(gsl_rng **r, int seed) 
{
	gsl_rng_env_setup();
	const gsl_rng_type *T;
	T = gsl_rng_mt19937;
	*r = gsl_rng_alloc(T);
	gsl_rng_set(*r, seed);

	for (int i = 0; i < 10000; i++)
		gsl_rng_uniform(*r);
}
