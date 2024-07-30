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
    
    // Constructor
    
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
    bool checkStatus();
};

inline void spinSite::Flip()
{
    spin *= -1;
}

inline void spinSite::AddToCluster()
{
    in_cluster = true;
}

inline void spinSite::Reset()
{
    in_cluster = false;
}

inline void spinSite::Triangulate(int i, int j)
{
    position[0] = i;
    position[1] = j;
}

inline int spinSite::getX()
{
    return position[0];
}

inline int spinSite::getY()
{
    return position[1];
}

inline int spinSite::getSpin()
{
    return spin;
}

inline bool spinSite::checkStatus() // change to getStatus();
{
    return in_cluster;
}