#include "parameters.h"
#include "spin_class.h"

void spinSite::AssignValues()
{
    int r;
    r = rand() % 2;
    if (r == 0) {
        spin = -1;
    } else {
        spin = 1;
    }
    in_cluster = false; // Start off "not in cluster"
}

void spinSite::Flip()
{
    spin *= -1;
}

void spinSite::AddToCluster()
{
    if (in_cluster == false)
    {
        in_cluster = true;
    }
}

void spinSite::Reset()
{
    in_cluster = false;
}

void spinSite::Triangulate(int i, int j)
{
    position[0] = i;
    position[1] = j;
}

int spinSite::getX()
{
    return position[0];
}

int spinSite::getY()
{
    return position[1];
}

int spinSite::getSpin()
{
    return spin;
}

int spinSite::checkStatus()
{
    return in_cluster;
}