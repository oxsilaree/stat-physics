#include <iostream>
#include <vector>
#include <bitset>
#include <math.h>
#include <memory>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <chrono>
#include <stack>
#include <random>
#include <list>
#include <vector>
#include "parameters.h"

using namespace std;

class spinSite
{// Access specifiers
private: // Data members
    
    int spin;
    bool in_cluster;
    int position[2];
    
public: // Member functions
    
    // Change data members
    void AssignValues();
    void Flip(); // changes spin
    void AddToCluster(); // changes in_cluster
    void Reset(); // changes in_cluster
    void Triangulate(int i, int j); // changes position

    // Get data members
    int getX();
    int getY();
    int getSpin();
    int checkStatus();
};