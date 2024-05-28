#include <iostream>
#include <stdio.h>
#include <fstream>
#include <math.h>
#include <time.h>
#include <chrono>
#include <stack>
#include <random>
#include <list>

using namespace std;

// Global variables

const int L = 64; // This is the only number we need to change
const int N = L*L;
int dim = 2;
const int nn_max = 6; // This means that we have both nearest neighbours and next-nearest neighbours in 1D

bool getlattice = true;
bool getdata = false;

const int blocks = 8;
const int sweeps = 128; //
const int steps = 16; // Average cluster size 

list<double> E;
list<double> M;
list<double> C;
list<double> X;

float J = 1;       // This is J1
double kappa;
double T;

string kappastr;
string Tstr;

double beta;
double padd1;
double padd2;

stack<int> step_counts;

class SpinCluster
{ // Here, let's declare a new class with a spin value and a criterion for whether it is already in the cluster or not

    // Access specifier
public:
    // Data members
    int spin;
    bool in_cluster;

    // Member functions
    void Flip()
    {
        spin *= -1;
    }
    void AddToCluster()
    {
        if (in_cluster == false)
        {
            in_cluster = true;
        }
    }
    void Reset()
    {
        in_cluster = false;
    }
};
SpinCluster lattice[N]; // Try signed char in the future
list<int> lattices;
int nnTable[N * nn_max]; // Make neighbor table

//FUNCTIONS
void Initialize(SpinCluster[], int);
double drandom();
int GetNeighbor(int, int);
void Step(int[]);
// void Annealing(int[]);
void ClusterAlgo(int[]);
double TotalEnergy(SpinCluster[], int[], int);
int TotalMag(SpinCluster[], int);
void LoadData();
void StepCounter();

int main()
{
  srand(time(NULL));
  cout << "input kappa:";
  cin >> kappa;
  cout << "input temperature:";
  cin >> T;

  beta = 1/T;
  cout << "starting now!!!\n";

  for (int i = 0; i < N; i++)
  {
    for (int j = 0; j < nn_max; j++)
    {
        nnTable[j * N + i] = GetNeighbor(i, j); 
        // cout << nnTable[j*L + i];
    }
  }

  Initialize(lattice, N); // Initialization


  for (int i = 0; i < N; i++) 
    {   // Store the lattice spin information (±1) to visualize later
        lattices.push_back(lattice[i].spin); 
    }

  ClusterAlgo(nnTable);
  
  cout << "so far so good!!!\n";
  LoadData();


  return 0;
}

void Initialize(SpinCluster lattice[], int N)
{
    int i, r;
    for (i = 0; i < N; i++)
    { // This fills the lattice
        r = rand() % 2;
        // cout << r;
        if (r == 0)
        {
            lattice[i].spin = -1;
        }
        else //if (r == 1)
        {
            lattice[i].spin = 1;
        }
        cout << lattice[i].spin << endl;
        lattice[i].in_cluster = false; // Initializes all the spins as NOT in cluster
    }
}

double drandom()
{
    return (double)rand()/RAND_MAX;
}

int GetNeighbor(int i, int pos) //i for seed index
{ // 0 -> left, 1 -> right, 2 -> up, 3 -> down, 4 -> up(NNN), 5 -> down(NNN)
    int l;
    if (pos == 0) {
        if (i % L == 0){
            l = i + (L-1);
        } else {
            l = i - 1;  
        }
    }
    if (pos == 1) {
        if (i % L == L - 1){
            l = i - (L-1);
        } else {
            l = i + 1;  
        }
    }
    if (pos == 2) { // up
        if (i < L){
            l = i + L*(L-1);
        } else {
            l = i - L;  
        }
    }
    if (pos == 3) { // down
        if (i >= L*(L-1)){
            l = i - L*(L-1);
        } else {
            l = i + L;  
        }
    }
    if (pos == 4) { // up NNN
        if (L <= i && i < 2*L){
            l = i + L*(L-2);
        } else {
            l = i - 2*L;  
        }
    }
    if (pos == 5) { // down NNN
        if (L*(L-2) <= i && i < L*(L-1)){
            l = i - L*(L-2);
        } else {
            l = i + 2*L;  
        }
    }
    return l;
}

void ClusterAlgo(int nnTable[]) 
{
    int i, m;
    double Mag, Ene, E1, E2, M1, M2, Mag_abs, M1_abs;
    padd1 = 1 - exp(-2 * beta * J);
    padd2 = 1 - exp(-2 * beta * J * kappa);
    
    for (i = 0; i < sweeps * 16; i++)
      {       // Burn in
        Step(nnTable);
      }

    for (m = 0; m < blocks; m++)
    {
        // measurements
        E1 = E2 = M1 = M2 = M1_abs = 0;
        for (i = 0; i < sweeps; i++)
        {
            for (int j = 0; j < steps; j++) {
                Step(nnTable); // Consider making one 'sweep' as a number of steps, where we choose it as after each spin has had one opportunity on average to flip
            }
            Ene = TotalEnergy(lattice, nnTable, L); // Look into changing when these are done for optimization
            Mag = TotalMag(lattice, L);
            Mag_abs = abs(Mag);
            M1_abs = M1_abs + Mag_abs;
            M1 = M1 + Mag;
            E1 = E1 + Ene;
            M2 = M2 + (Mag * Mag);
            
            E2 = E2 + (Ene * Ene);
        }

        E.push_back(E1 / (sweeps)); // Total Energy of lattice
        M.push_back(M1 / (sweeps)); // Total Magnetization of lattice
        C.push_back(((E2 / sweeps) - pow((E1 / (sweeps)),2)) * pow(beta,2)); // Specific heat  (not per spin) C = var(E)/(k_B T)
        X.push_back(((M2 / sweeps) - pow((M1 / (sweeps)),2)) * pow(beta,2)); // Susceptibility (not per spin)

        for (int i = 0; i < N; i++) 
            {   // Store the lattice spin information (±1) to visualize later
                lattices.push_back(lattice[i].spin); 
            }
    }
}

void Step(int nnTable[])
{

    int i, k, sp, oldspin, newspin, current, nn;
    int stacker[N];

    i = rand() % N;

    stack<int> cluster;

    stacker[0] = i;

    sp = 1;
    
    lattice[i].AddToCluster();
    cluster.push(i);
    int step_count = 0;
    while (sp)
    {

        /* Pull a site off the stack from the nearest neighbours*/

        current = stacker[--sp];
        oldspin = lattice[current].spin;
        newspin = -oldspin;

        /* Check the neighbours using neighbour table */

        for (int j = 0; j < nn_max; j++)
        { 
            nn = nnTable[N * j + current];
            if (j <= 3) // Nearest neighbours
            {
                
                if (lattice[nn].spin == oldspin && lattice[nn].in_cluster == false)
                {
                    if (drandom() <= padd1)
                    {
                        stacker[sp++] = nn;
                        lattice[nn].AddToCluster();
                        cluster.push(nn);
                        step_count += 1;
                    }
                }
            }
            else if (j >= 4)
            { // Next nearest neighbours
                if (lattice[nn].spin == newspin && lattice[nn].in_cluster == false)
                {
                    if (drandom() <= padd2)
                    {
                        stacker[sp++] = nn;
                        lattice[nn].AddToCluster();
                        cluster.push(nn);
                        step_count += 1;
                    }
                }
            }
        }

        step_counts.push(step_count);
    }
    /* Go over all the spins in the stack and flip if they are in the cluster */
    while (!cluster.empty()) {
        k = cluster.top();
        lattice[k].Flip();
        lattice[k].Reset();
        cluster.pop();
    }
}

double TotalEnergy(SpinCluster lattice[], int nnTable[], int N)
{
    double energy = 0;
    int left, right, up, down;   // nearest neighbours
    int up2, down2; // next-nearest neighbours
    int J1 = 1;
    for (int i = 0; i < N; i++)
    {
        left = nnTable[N * 0 + i];
        right = nnTable[N * 1 + i];
        up = nnTable[N * 2 + i];
        down = nnTable[N * 3 + i];
        up2 = nnTable[N * 4 + i];
        down2 = nnTable[N * 5 + i];
        

        // Assuming J = 1
        // energy += lattice[i].spin * -(left + right);
        // energy += lattice[i].spin * kappa * (left2 + right2);
        energy += lattice[i].spin * (-1*J1) * ((lattice[left].spin + lattice[right].spin + lattice[up].spin + lattice[down].spin) - (kappa * (lattice[up2].spin + lattice[down2].spin)));
    }
    energy *= 0.5;
    return energy;
}

int TotalMag(SpinCluster lattice[], int N)
{
    int i;
    int mag = 0;
    for (i = 0; i < N; i++)
        mag += lattice[i].spin;
    return mag;
}


void LoadData() 
{
    if (getlattice == true) {
        ofstream lattice_data;
        string directoryPath = "/Users/shanekeiser/Documents/Spring 2024/Machta Spring/2D_ANNNI/2d_lattices/";
        lattice_data.open(directoryPath + "2d_ANNNI_annealing_visual.csv");
        
        for (list<int>::iterator it=lattices.begin(); it != lattices.end(); ++it)
        {
            lattice_data << *it << ",";
        }
    lattice_data.close();
    }

    if (getdata == true){
        ofstream mdata; // make a file to write to
        string directoryPath = "/Users/shanekeiser/Documents/Spring 2024/Machta Spring/2D_ANNNI/2d_data/";
        mdata.open(directoryPath + "1d_ANNNI_annealing.csv");
        mdata << "E,";
        
        for (list<double>::iterator it=E.begin(); it != E.end(); ++it)
        {
            mdata << *it << ",";
        }
        mdata << "\n";
        mdata << "M,";

        for (list<double>::iterator it=M.begin(); it != M.end(); ++it)
        {
            mdata << *it << ",";
        }
        mdata << "\n";
        mdata << "C,";

        for (list<double>::iterator it=C.begin(); it != C.end(); ++it)
        {
            mdata << *it << ",";
        }
        mdata << "\n";
        mdata << "X,";

        for (list<double>::iterator it=X.begin(); it != X.end(); ++it)
        {
            mdata << *it << ",";
        }
        mdata << "\n";
        
    mdata.close();
    }
}

void StepCounter() 
{
    float avg_steps; // Number of spins in cluster at each step
    int len_step_counts = step_counts.size();
    while (!step_counts.empty())
    {
        avg_steps += step_counts.top();
        step_counts.pop();
    }
    avg_steps /= len_step_counts;
    std::cout << avg_steps;
}