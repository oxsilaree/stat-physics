/*
Most parameters are self-explanatory. Array_index is mostly obsolete, but is part of
creating a 3D neighbor table in case you don't have an external table.
*/

#define LENGTH 10
#define LENGTH_3 (LENGTH *LENGTH *LENGTH)
#define array_index(i, j, k) (k + j *LENGTH + i *LENGTH *LENGTH)
#define CULLING_FRAC 0.15
#define MAX_BETA 5
#define PI 3.14159265358979323846
#define NUM_THREADS 40
