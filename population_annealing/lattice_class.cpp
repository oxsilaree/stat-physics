#include "lattice_class.h"

/*
void Lattice::initializeSites()
{
    int i, j;
    for (i = 0; i < LEN; i++) 
    {
        for (j = 0; j < LEN; j++) 
        {
            lattice_object[i][j].AssignValues();
        }
    }
}
*/

// Constructor definition
Lattice::Lattice()
    {
    // Initialize other members if needed
    energy = 0;
    mag = 0;
    abs_mag = 0;
    avg_cluster_size = 0;
    spec_heat = 0.0;
    suscep = 0.0;
    lattice_number = 0;
    lattice_object = vector<vector<spinSite> >(LEN, vector<spinSite>(LEN));
    }  


void Lattice::initializeSites()
{

    for (int i = 0; i < LEN; i++)
    {
        for (int j = 0; j < LEN; j++)
        {
            // cout << " ---i = " << i << ", j = " << j <<"--- \n";
            spinSite* starter = getSpinSite(i,j);
            starter->AssignValues();
            // cout << "Spin is " << starter->getSpin() << ". \n";
        }
    }

}

spinSite* Lattice::getSpinSite(int row, int col)
{
    // cout <<"---row = " << row << ", col = " << col << "---\n";
    return &lattice_object[(int)row][(int)col];
}

void Lattice::doBurnIn(int neighbor_table[LEN][LEN][nn_max][dim], double T)
{
    double Beta, padd1, padd2;
    Beta = 1/T;
    padd1 = 1 - exp(-2 * Beta * J);
    padd2 = 1 - exp(-2 * Beta * J * kappa);
    for (int i = 0; i < SWEEPS*STEPS*4; i++)
    {
        doBurnInStep(neighbor_table, padd1, padd2);
    }
}

void Lattice::doBurnInStep(int neighbor_table[LEN][LEN][nn_max][dim], double padd1, double padd2)
{ // Basically doStep without the wrapping checks
    int i, j,   lx, ly,     oldspin, newspin,   current_x, current_y,   nn_i, nn_j; // Lattice indices
    int root_x, root_y,     coord_x, coord_y,   pos_update_x, pos_update_y; // Coordinates for wrapping criteria
    int old_x, old_y,       new_x, new_y;       // More coordinates for wrapping criteria 
    int wrapcounter, sp, cluster_size; // Counters
    double randnum; // For checking whether to add to cluster
    bool wrapping_crit = false;
    // spinSite root_spin, current_spin, neighbor_spin;
    vector<pair<int, int> > stacker;
    // stack<int> cluster_x;
    // stack<int> cluster_y;

    // Pick a random spin (x then y coord) and make it the origin (i.e. "root/seed/(0,0)")
    i = rand() % LEN;
    j = rand() % LEN;
    // cout << "--- i = " << i << ", j = " << j << " ---\n";
    spinSite* root_spin = getSpinSite(i,j);
    // cout << root_spin->getSpin(); ///////// DEBUGGING PRINTER
    root_spin->AddToCluster();
    root_spin->Triangulate(0,0);

    // cluster_x.push(i);
    // cluster_y.push(j);
    stacker.push_back(make_pair(i, j));
    cluster_size = 1;

    sp = 1;
    wrapcounter = 0;
    while (sp)
    {   /* Pull a site off the stack from the nearest neighbours*/
        --sp;
        current_x = stacker.back().first;
        current_y = stacker.back().second;
        stacker.pop_back();
        spinSite* current_spin = getSpinSite(current_x, current_y);
        oldspin = current_spin->getSpin();
        newspin = -oldspin;
        coord_x = current_spin->getX();
        coord_y = current_spin->getY();

        /* Check the neighbours using neighbour table */
        for (int k = 0; k < nn_max; k++)
        { 
            nn_i = neighbor_table[current_x][current_y][k][0];
            nn_j = neighbor_table[current_x][current_y][k][1];
            spinSite* neighbor_spin = getSpinSite(nn_i, nn_j);
            if (neighbor_spin->checkStatus() == false) // Nearest neighbours
            {
                randnum = (double)rand()/RAND_MAX;   
                if (neighbor_spin->getSpin() == oldspin) 
                {
                                                                                                 
                    if (randnum <= padd1)
                    {
                        if (k == 0){
                            pos_update_x = coord_x - 1;
                            pos_update_y = coord_y;
                        } else if (k == 1) {
                            pos_update_x = coord_x + 1;
                            pos_update_y = coord_y;
                        } else if (k == 2) {
                            pos_update_x = coord_x;
                            pos_update_y = coord_y + 1;
                        } else if (k == 3) {
                            pos_update_x = coord_x;
                            pos_update_y = coord_y - 1;
                        }
                        sp += 1;
                        stacker.push_back(make_pair(nn_i, nn_j));
                        neighbor_spin->AddToCluster();
                        neighbor_spin->Triangulate(pos_update_x, pos_update_y); 
                        // cluster_x.push(nn_i);
                        // cluster_y.push(nn_j);
                        cluster_size += 1;
                    }
                }
                // Next nearest neighbours
                if (neighbor_spin->getSpin() == newspin)
                {
                    if (randnum <= padd2)
                    {
                        if (k == 4) {
                            pos_update_x = coord_x;
                            pos_update_y = coord_y + 2;
                        } else if (k == 5) {
                            pos_update_x = coord_x;
                            pos_update_y = coord_y - 2;
                        }
                        sp += 1;
                        stacker.push_back(make_pair(nn_i, nn_j));
                        neighbor_spin->AddToCluster();
                        neighbor_spin->Triangulate(pos_update_x, pos_update_y); 
                        // cluster_x.push(nn_i);
                        // cluster_y.push(nn_j);
                        cluster_size += 1;
                    }
                }
            }
        }
    }
    /* Go over all the spins in the stack and flip if they are in the cluster */
    while (!stacker.empty()) { //!cluster_x.empty()
        //lx = cluster_x.top();
        //ly = cluster_y.top();
        lx = stacker.back().first;
        ly = stacker.back().second;
        spinSite* flipper = getSpinSite(lx,ly);
        flipper->Flip();
        flipper->Reset();
        stacker.pop_back();
        //cluster_x.pop();
        //cluster_y.pop();
    }
}



void Lattice::doStep(int neighbor_table[LEN][LEN][nn_max][dim], double padd1, double padd2)
{
    int i, j,   lx, ly,     oldspin, newspin,   current_x, current_y,   nn_i, nn_j; // Lattice indices
    int root_x, root_y,     coord_x, coord_y,   pos_update_x, pos_update_y; // Coordinates for wrapping criteria
    int old_x, old_y,       new_x, new_y;       // More coordinates for wrapping criteria 
    int wrapcounter, sp, cluster_size; // Counters
    double randnum; // Number to decide if add to cluster
    bool wrapping_crit = false;
    vector<pair<int, int> > stacker;
    stack<int> cluster_x;
    stack<int> cluster_y;

    // Pick a random spin (x then y coord) and make it the origin (i.e. "root/seed/(0,0)")
    i = rand() % LEN;
    j = rand() % LEN;
    spinSite* root_spin = getSpinSite(i,j);
    // cout << root_spin->getSpin(); ///////// DEBUGGING PRINTER
    root_spin->AddToCluster();
    root_spin->Triangulate(0,0);
    // lattice_object[i][j].AddToCluster();
    // lattice_object[i][j].Triangulate(0,0);
    cluster_x.push(i);
    cluster_y.push(j);
    stacker.push_back(make_pair(i, j));
    cluster_size = 1;

    sp = 1;
    wrapcounter = 0;
    while (sp)
    {   /* Pull a site off the stack from the nearest neighbours*/
        --sp;
        current_x = stacker.back().first;
        current_y = stacker.back().second;
        stacker.pop_back();
        spinSite* current_spin = getSpinSite(current_x, current_y);
        oldspin = current_spin->getSpin();
        newspin = -oldspin;
        coord_x = current_spin->getX();
        coord_y = current_spin->getY();
        //oldspin = lattice_object[current_x][current_y].getSpin();
        //newspin = -oldspin;
        //coord_x = lattice_object[current_x][current_y].getX();
        //coord_y = lattice_object[current_x][current_y].getY();

        /* Check the neighbours using neighbour table */

        for (int k = 0; k < nn_max; k++)
        {   
            nn_i = neighbor_table[current_x][current_y][k][0];
            nn_j = neighbor_table[current_x][current_y][k][1];
            if (nn_i < 0 || nn_j < 0 || nn_i >= LEN || nn_j >= LEN)
            {
                cout << "---nn_i = " << nn_i << ", nn_j = " << nn_j << "---\n";
            }
            spinSite* neighbor_spin = getSpinSite(nn_i,nn_j);
            if (neighbor_spin->checkStatus() == true && wrapping_crit == false) {

                old_x = neighbor_spin->getX();
                old_y = neighbor_spin->getY();
                if (k == 0){ // left bond
                            pos_update_x = coord_x - 1;
                            pos_update_y = coord_y;
                        } else if (k == 1) { // right bond
                            pos_update_x = coord_x + 1;
                            pos_update_y = coord_y;
                        } else if (k == 2) { // up bond
                            pos_update_x = coord_x;
                            pos_update_y = coord_y + 1;
                        } else if (k == 3) { // down bond
                            pos_update_x = coord_x;
                            pos_update_y = coord_y - 1;
                        } else if (k == 4) { // up2 bond
                            pos_update_x = coord_x;
                            pos_update_y = coord_y + 2;
                        } else if (k == 5) { // down2 bond
                            pos_update_x = coord_x;
                            pos_update_y = coord_y - 2;
                        }
                neighbor_spin->Triangulate(pos_update_x, pos_update_y);
                new_x = neighbor_spin->getX();
                new_y = neighbor_spin->getY();

                if (old_x != new_x || old_y != new_y) { // Check if either of the new-coordinates (rel. to root spin) is different from previous inclusion.
                    wrapping_crit = true; // If it is different, then we have 'wrapped' around the lattice.    
                    wrapcounter++;
                    // cout << "Wrapping has occurred.\n";
                    // wrapcounts.push(wrapcounter);                   
                }
            }
            
            else if (neighbor_spin->checkStatus() == false) // Nearest neighbours
            {
                randnum = (double)rand()/RAND_MAX;
                if (neighbor_spin->getSpin() == oldspin) 
                {                                                                                      
                    if (randnum <= padd1)
                    {
                        if (k == 0){
                            pos_update_x = coord_x - 1;
                            pos_update_y = coord_y;
                        } else if (k == 1) {
                            pos_update_x = coord_x + 1;
                            pos_update_y = coord_y;  
                        } else if (k == 2) {
                            pos_update_x = coord_x;
                            pos_update_y = coord_y + 1;
                        } else if (k == 3) {
                            pos_update_x = coord_x;
                            pos_update_y = coord_y - 1;
                        }
                        sp += 1;
                        stacker.push_back(make_pair(nn_i, nn_j));
                        neighbor_spin->AddToCluster();
                        neighbor_spin->Triangulate(pos_update_x, pos_update_y); 
                        cluster_x.push(nn_i);
                        cluster_y.push(nn_j);
                        cluster_size += 1;
                    }
                }
            
                // Next nearest neighbours
                if (neighbor_spin->getSpin() == newspin)
                {   
                    if (randnum <= padd2)
                    {
                        if (k == 4) {
                            pos_update_x = coord_x;
                            pos_update_y = coord_y + 2;
                        } else if (k == 5) {
                            pos_update_x = coord_x;
                            pos_update_y = coord_y - 2;
                        }
                        sp += 1;
                        stacker.push_back(make_pair(nn_i, nn_j));
                        neighbor_spin->AddToCluster();
                        neighbor_spin->Triangulate(pos_update_x, pos_update_y); 
                        cluster_x.push(nn_i);
                        cluster_y.push(nn_j);
                        cluster_size += 1;
                    }
                }
            }
        }
    }
    /* FIGURE THIS OUT LATER
    if (wrapping_crit == true){
        cluster_sizes.push_back(-1*cluster_size); // We use negative cluster size to indicate clusters that have wrapped,
    } else {                                 // As a slick way to keep everything in one array(stack)
        cluster_sizes.push_back(cluster_size);
    }
    */
    /* Go over all the spins in the stack and flip if they are in the cluster */
    while (!cluster_x.empty()) {
        lx = cluster_x.top();
        ly = cluster_y.top();
        lattice_object[lx][ly].Flip();
        lattice_object[lx][ly].Reset();
        cluster_x.pop();
        cluster_y.pop();
    }
}

void Lattice::doSweep(int neighbor_table[LEN][LEN][nn_max][dim], double T)
{
    int Ene, E1, E2,    Mag, M1, M2, M1_abs, Mag_abs, Mag_sq;
    double Beta, padd1, padd2;
    Beta = 1/T;
    padd1 = 1 - exp(-2 * Beta * J);
    padd2 = 1 - exp(-2 * Beta * J * kappa);
    
    E1 = E2 = M1 = M2 = M1_abs = 0;
        for (int j = 0; j < STEPS; j++) {
            doStep(neighbor_table, padd1, padd2); // Consider making one 'sweep' as a number of steps, where we choose it as after each spin has had one opportunity on average to flip
        }
        updateTotalEnergy(neighbor_table); // Take data after each sweep
        updateTotalMag();
        Ene = getTotalEnergy();
        Mag_abs = abs(Mag);
        Mag_sq = pow(Mag,2);

        M1_abs = M1_abs + Mag_abs;
        M1 = M1 + Mag;
        E1 = E1 + Ene;
        M2 = M2 + Mag_sq;
        E2 = E2 + (Ene * Ene);
}

void Lattice::doWolffAlgo(int neighbor_table[LEN][LEN][nn_max][dim], double T, int num_sweeps)
{
    // Initialize and Burn-in
    doBurnIn(neighbor_table, T);
    for (int m = 0; m < BLOCKS; m++)
    {
        for (int i = 0; i < num_sweeps; i++)
        {
            doSweep(neighbor_table, T);
        }
    }
}

void Lattice::updateTotalEnergy(int neighbor_table[LEN][LEN][nn_max][dim])
{
    int left_x, right_x, up_x, down_x;   // nearest neighbours
    int left_y, right_y, up_y, down_y;
    int up2_x, down2_x; // next-nearest neighbours
    int up2_y, down2_y;

    
    for (int i = 0; i < LEN; i++)
    {
        for (int j = 0; j < LEN; j++)
        {
            left_x = neighbor_table[i][j][0][0];
            left_y = neighbor_table[i][j][0][1];
            right_x = neighbor_table[i][j][1][0];
            right_y = neighbor_table[i][j][1][1];
            up_x = neighbor_table[i][j][2][0];
            up_y = neighbor_table[i][j][2][1];
            down_x = neighbor_table[i][j][3][0];
            down_y = neighbor_table[i][j][3][1];
            up2_x = neighbor_table[i][j][4][0];
            up2_y = neighbor_table[i][j][4][1];
            down2_x = neighbor_table[i][j][5][0];
            down2_y = neighbor_table[i][j][5][1];

            spinSite* mid_s = getSpinSite(i,j);
            spinSite* left_s = getSpinSite(left_x,left_y);
            spinSite* right_s = getSpinSite(right_x,right_y);
            spinSite* up_s = getSpinSite(up_x, up_y);
            spinSite* down_s = getSpinSite(down_x, down_y);
            spinSite* up2_s = getSpinSite(up2_x, up2_y);
            spinSite* down2_s = getSpinSite(down2_x, down2_y);            


        energy += (mid_s->getSpin())*(-1*J) * (((left_s->getSpin())+(right_s->getSpin()) \
                                            + (up_s->getSpin()) + (down_s->getSpin())) \
                                            - (kappa * ((up2_s->getSpin()) + (down2_s->getSpin()))));
        /*
        energy += lattice_object[i][j].getSpin() * (-1*J) \
                * ((lattice_object[left_x][left_y].getSpin() + lattice_object[right_x][right_y].getSpin() \
                + lattice_object[up_x][up_y].getSpin() + lattice_object[down_x][down_y].getSpin()) \
                - (kappa * (lattice_object[up2_x][up2_y].getSpin() + lattice_object[down2_x][down2_y].getSpin())));
        */
        }  
    }   
    energy *= 0.5;
}

void Lattice::updateTotalMag()
{
    int i, j;
    mag = 0;
    for (i = 0; i < LEN; i++){
        for (j = 0; j < LEN; j++){
            mag += lattice_object[i][j].getSpin();
        }
    }
}

int Lattice::getTotalEnergy()
{
    return energy;
}

int Lattice::getTotalMag()
{
    return mag;
}