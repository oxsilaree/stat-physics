#define kappa 2        // -J2/J1, ratio of coupling constants. kappa = 0 --> regular Ising
#define J 1
#define nn_max 6        // No. of neighbors
#define dim 2           // Dimensions
#define L 16            // System side length
#define N (pow(L,dim))  // System size (L^dim)
#define T_init 5        // In pop. annealing we will lower this as we run the simulation
#define PI 3.14159265358979323846
#define NUM_THREADS 40

#define blocks 10
#define sweeps 20
#define steps 50