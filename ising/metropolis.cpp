#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>

using namespace std;


// Declare functions
void metropolis();
int deltaU( int );
void initialize( int [],         int );
void sweep(      int [], double, int );
int total_energy(int [],         int );
int total_mag( int [],           int );

void wolffcluster(int);
void wolff();
void clusterAdder(int);
    

// Declare global variables

const int L = 45; // This is the only number we need to change
const int N = L*L;
int lattice[N]; // Try signed char in the future
int buffer[N]; // for Wolff subroutine

bool getlattice = 1;
bool getdata = 1;

const double T = 2.27;
const double beta = 1/T;

double wolffprob = 1 - exp(-2*beta); // J = 1

const int blocks = 16;
const int sweeps = N * 8;

double E[blocks];
double M[blocks];
double C[blocks];
double X[blocks];

int lattices[N*blocks];

int main() {

    // RUN METROPOLIS ALGORITHM
    metropolis();

    // WRITE LATTICE FILE
    if (getlattice == 1) {
        ofstream lattice_data;
        lattice_data.open("lattice_visual_" + to_string(L) + "L.csv");
        
        for (int k = 0; k < blocks*N; k++) {
        lattice_data << lattices[k] << ',';
            }
        lattice_data.close();
    }
    
    // WRITE DATA FILE FOR E,M,C,X

    if (getdata == 1){
        ofstream mdata; // make a file to write to
        mdata.open("new_data_" + to_string(L) + "L.csv");
        mdata << "E, M, C, X" << "\n";
        for (int i = 0; i < blocks; i++) {
            mdata << E[i] << ", ";
            mdata << M[i] << ", ";
            mdata << C[i] << ", ";
            mdata << X[i] << "\n";
        }
    mdata.close();

    
    }
    return 0;
}

void wolff() {
    srand(time(NULL));
    initialize(lattice, N);

    
    int index = rand() % N;
    wolffcluster(index);

    for (int i = 0; i < N; i++) {
        cout << buffer[i];
    }
}

void wolffcluster(int index) {
    
    for (int i = 0; i < N; i++) {
        buffer[i] = 0;
    }

    buffer[index] = -lattice[index];
    
    clusterAdder(index);

}

void clusterAdder(int i) {
    bool checkDone = false;
    int left,right,top,bottom;
    while (checkDone == false) {

        if (i % L == 0){ // Left BC
            left = i + (L-1); 
            } else { left = i - 1; }

        if (i % L == L - 1){ // Right BC
            right = i - (L-1); 
            } else { right = i + 1; }

        if (i <= L - 1) { // Top BC
            top = i + (N - L); 
            } else { top = i - L; }

        if (i >= N - L) { // Bottom BC
            bottom = i - (N-L);
        } else { bottom = i + L; }
        

        int surrounds[4] = { left, right, top, bottom };
        for (int j = 0; j < 4; j++) {
            if (lattice[surrounds[j]] == lattice[i] && buffer[surrounds[j]] == 0 && (rand()%1000000)/1000000.0 <= wolffprob) {
                buffer[surrounds[j]] = -lattice[surrounds[j]];
                clusterAdder(surrounds[j]);
            }
            
        }
    }
}

void metropolis() {
    srand(time(NULL));
    
    int i, m;
    double Mag, Ene, E1, E2, M1, M2;

    // initialize the SYSTEM
    initialize(lattice, N);		

    // equilibrating the SYSTEM 
    for (i = 0 ; i < sweeps * 4 ; i++){ // Burn in
        sweep(lattice, beta, N);
        }
    
    for (m = 0; m < blocks; m++) {
        // measurements
        E1 = E2 = M1 = M2 = 0;
        for (i = 0 ; i < sweeps ; i++){
            sweep(lattice, beta, N);
                
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
// FUNCTIONS USED.
/* 1. Initialize the lattice on the Lattice with periodic boundary conditions */
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


/* 2. Monte Carlo Moves with the Metropolis Algorithm */
	void sweep( int lattice[], double beta, int N ){
	int i, ipick, Ediff;
  	double cost;
  		for (i = 0 ; i < N ; i++){
		ipick =  rand() % N;	  
     	Ediff = deltaU(ipick);
        cost = exp(-Ediff * beta);
			  
			if (Ediff < 0) {
                lattice[ipick] = -lattice[ipick];
            } else {
			  	 if ((rand()%1000000)/1000000.0 <= cost) {
                    lattice[ipick]= -lattice[ipick];
                    } else {
                        lattice[ipick] = lattice[ipick];
                    }
  			}
	    }
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
