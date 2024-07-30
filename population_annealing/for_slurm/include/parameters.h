// #define kappa 0.0
#define J 1.0
#define NN_MAX 6            // No. of neighbors
#define DIM 2               // Dimensions
#define LEN 16           // System side length
// #define N (pow(L,dim))  // System size (L^dim) -- probably not needed. The gsl_rng library uses N, oops
// #define T_INIT 3.0        // In pop. annealing we will lower this as we run the simulation
// #define T_FINAL 1.0    // Set minimum T to 0.25 so our beta doesn't blow up so much
#define T_ITER 100.0
#define PI 3.14159265358979323846
#define NUM_THREADS 500
#define INIT_POP_SIZE 500
#define CULLING_FRAC 0.02
#define MAX_BETA 2.5

#define BLOCKS 1 // This number is redundant now.
#define SWEEPS 4
#define STEPS 16
