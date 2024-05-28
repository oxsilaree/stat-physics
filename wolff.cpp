#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <chrono>
using namespace std;
using namespace std::chrono;

// Declare global variables

const int L = 50; // This is the only number we need to change
const int N = L*L;
int lattice[N]; // Try signed char in the future

bool getlattice = 1;
bool getdata = 1;

const double T = 2.27;
const double beta = 1/T;

const int blocks = 16;
const int sweeps = N * 8;

double E[blocks];
double M[blocks];
double C[blocks];
double X[blocks];

int lattices[N*blocks];

void wolff();
void initialize(int[], int);
void step();
double drandom();
int total_energy(int [],         int );
int total_mag( int [],           int );
int deltaU( int );


int main() {
    auto start = high_resolution_clock::now();

    wolff();

    auto stop = high_resolution_clock::now();

    // WRITE DATA FILE FOR E,M,C,X

    if (getdata == 1){
        string directoryPath = "/Users/shanekeiser/Documents/researchMachta/code/spring/data/";
        ofstream mdata; // make a file to write to
        mdata.open("wolff_data" + to_string(L) + "L.csv");
        mdata << "E, M, C, X" << "\n";
        for (int i = 0; i < blocks; i++) {
            mdata << E[i] << ", ";
            mdata << M[i] << ", ";
            mdata << C[i] << ", ";
            mdata << X[i] << "\n";
        }
    mdata.close();
    
    }
    
    auto duration = duration_cast<microseconds>(stop - start);
    cout << duration.count() << endl;
    return 0;
}


void wolff() {

    srand(time(NULL));
    
    int i, m;
    double Mag, Ene, E1, E2, M1, M2;

    #define XNN 1
    #define YNN L


    // initialize the SYSTEM
    initialize(lattice, N);	
    

    // equilibrating the SYSTEM 
    for (i = 0 ; i < sweeps * 4 ; i++){ // Burn in
        step();
        }
    
    for (m = 0; m < blocks; m++) {
        // measurements
        E1 = E2 = M1 = M2 = 0;
        for (i = 0 ; i < sweeps ; i++){
            step();
                
            Mag  = total_mag(lattice, N);
            Ene  = total_energy(lattice, N);
            
            M1 = M1   + Mag;
            E1 = E1   + Ene;
            M2 = M2   + Mag*Mag ;
            E2 = E2   + Ene*Ene;
            }

        E[m]  = E1/(N*sweeps);
        M[m]  = M1/(N*sweeps);
        C[m]  = ( E2/sweeps - E1*E1/(sweeps*sweeps) )*beta*beta/(N);
        X[m]  = ( M2/sweeps - M1*M1/(sweeps*sweeps) )*beta*beta/(N);
        for (int n = 0; n < N; n++) {
            lattices[m*N + n] = lattice[n];
        }
    }


}

void initialize( int lattice[], int N ){
srand(time(0)); // Random seed
    int i,r;
    for (i = 0; i < N; i++) { // This fills the lattice
        r = rand() % 2;
        // cout << r;
        if (r == 0) {
            lattice[i] = -1;
        } else if (r == 1) { 
            lattice[i] = 1;
        }
    }
}


void step() {

    int i;
    int sp;
    int oldspin,newspin;
    int current, nn;
    int stack[N];

    int J = 1;
    double padd = 1 - exp(-2*beta*J);

    // i = N*(double)rand()/RAND_MAX;
    i = rand() % N;


    stack[0] = i;
    sp = 1;
    oldspin = lattice[i];
    newspin = -lattice[i];
    lattice[i] = newspin;

    while (sp) {

        /* Pull a site off the stack */

        current = stack[--sp];

        /* Check the neighbours */

        if ((nn=current+XNN)>=N) nn -= N;
        if (lattice[nn]==oldspin)
            if (drandom()<padd) {
            stack[sp++] = nn;
            lattice[nn] = newspin;
            }
        if ((nn=current-XNN)<0) nn += N;
        if (lattice[nn]==oldspin)
            if (drandom()<padd) {
            stack[sp++] = nn;
            lattice[nn] = newspin;
            }
        if ((nn=current+YNN)>=N) nn -= N;
        if (lattice[nn]==oldspin)
            if (drandom()<padd) {
            stack [sp++] = nn;
            lattice[nn] = newspin;
            }
        if ((nn=current-YNN)<0) nn += N;
        if (lattice[nn]==oldspin)
            if (drandom()<padd) {
            stack[sp++] = nn;
            lattice[nn] = newspin;
            }
    }
}

double drandom() {
    return (double) rand()/RAND_MAX;
}

/* 3. Total Energy */
int total_energy( int lattice[], int N ){
    int energy = 0;
    int left,right,top,bottom; 
    for (int i = 0; i < N; i++) {
        if (i % L == 0){ // Left BC
            left = lattice[i + (L-1)]; 
            } else { left = lattice[i - 1]; }

        if (i % L == L - 1){ // Right BC
            right = lattice[i - (L-1)]; 
            } else { right = lattice[i + 1]; }

        if (i <= L - 1) { // Top BC
            top = lattice[i + (N - L)]; 
            } else { top = lattice[i - L]; }

        if (i >= N - L) { // Bottom BC
            bottom = lattice[i - (N-L)];
        } else { bottom = lattice[i + L]; }

        // Assuming J = 1
        energy += lattice[i] * -(top+bottom+left+right);
        
    }
    energy *= 0.5;
    return energy;
}


/* 4. Total Magnetization */
int total_mag( int lattice[], int N ){
	int i;
	double mag= 0;
        for (i = 0 ; i < N; i++) mag = mag + lattice[i];
        return(mag);
}

int deltaU(int i) { // Periodic boundary conditions
    int top, bottom, left, right;
    int E_diff;
    
    if (i % L == 0){ // Left BC
        left = lattice[i + (L-1)]; 
        } else { left = lattice[i - 1]; }

    if (i % L == L - 1){ // Right BC
        right = lattice[i - (L-1)]; 
        } else { right = lattice[i + 1]; }

    if (i <= L - 1) { // Top BC
        top = lattice[i + (N - L)]; 
        } else { top = lattice[i - L]; }

    if (i >= N - L) { // Bottom BC
        bottom = lattice[i - (N - L)];
    } else { bottom = lattice[i + L]; }

    // Assuming J = 1
    E_diff = 2*lattice[i]*(top+bottom+left+right);

    return E_diff;
}
