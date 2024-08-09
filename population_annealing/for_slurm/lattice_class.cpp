#include "lattice_class.h"

// Constructor definition
Lattice::Lattice(double kappa)
    {
    // Initialize other members if needed
    Lattice::kappa = kappa;
    Lattice::energy = 0;
    Lattice::mag = 0;
    Lattice::abs_mag = 0;
    Lattice::avg_cluster_size = 0;
    Lattice::avg_nowrap_cluster_size = 0;
    Lattice::wrap_counter = 0;
    Lattice::nowrap_counter = 0;
    Lattice::spec_heat = 0.0;
    Lattice::suscep = 0.0;
    Lattice::lattice_number = 0;
    Lattice::lattice_object = vector<vector<spinSite> >(LEN, vector<spinSite>(LEN));
    }


void Lattice::initializeSites(double *Beta)
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
    // Initialize and Burn-in
    if (*Beta != 0) { // 23/7/24 if we start at inf temp, no need to burn in. Just resample
        doBurnIn(*Beta);
    }
     
}



void Lattice::doBurnIn(double Beta)
{
    double padd1, padd2;
    padd1 = 1 - exp(-2 * Beta * J);
    padd2 = 1 - exp(-2 * Beta * J * kappa);
    for (int i = 0; i < SWEEPS*STEPS*4; i++)
    {
        doBurnInStep(&padd1, &padd2);
    }
}

void Lattice::doBurnInStep(double *padd1, double *padd2)
{ // Basically doStep without the wrapping checks
    int i, j,   lx, ly,     oldspin, newspin,   current_x, current_y,   nn_i, nn_j; // Lattice indices
    int root_x, root_y,     coord_x, coord_y,   pos_update_x, pos_update_y; // Coordinates for wrapping criteria
    int sp; // Counters
    double randnum; // For checking whether to add to cluster
    // spinSite root_spin, current_spin, neighbor_spin;
    vector<pair<int, int> > stacker;
    stacker.reserve(LEN*LEN);
    stack<int> cluster_x;
    stack<int> cluster_y;
    i = rand() % LEN;
    j = rand() % LEN;
    spinSite* root_spin = getSpinSite(i,j);
    root_spin->AddToCluster();
    root_spin->Triangulate(0,0);
    cluster_x.push(i);
    cluster_y.push(j);
    stacker.push_back(make_pair(i, j));

    sp = 1;
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
        for (int k = 0; k < NN_MAX; k++)
        { 
            nn_i = neighbor_table[current_x][current_y][k][0];
            nn_j = neighbor_table[current_x][current_y][k][1];
            spinSite* neighbor_spin = getSpinSite(nn_i, nn_j);
            bool neighbor_checked = neighbor_spin->checkStatus();
            if (neighbor_checked == false)  {
                randnum = static_cast<double>(rand()) / RAND_MAX;
                if ((neighbor_spin->getSpin() == oldspin && k <= 3 && randnum <= *padd1) ||
                    (neighbor_spin->getSpin() == newspin && k >= 4 && randnum <= *padd2)) 
                {
                    
                    int pos_update_x = coord_x, pos_update_y = coord_y;
                    if (k == 0) pos_update_x -= 1;
                    else if (k == 1) pos_update_x += 1;
                    else if (k == 2) pos_update_y += 1;
                    else if (k == 3) pos_update_y -= 1;
                    else if (k == 4) pos_update_y += 2;
                    else if (k == 5) pos_update_y -= 2;

                    sp += 1;
                    stacker.push_back(make_pair(nn_i, nn_j));
                    neighbor_spin->AddToCluster();
                    neighbor_spin->Triangulate(pos_update_x, pos_update_y);
                    cluster_x.push(nn_i);
                    cluster_y.push(nn_j);
                }
            }
        }
    }
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



void Lattice::doStep(double *padd1, double *padd2)
{
    int i, j,   lx, ly,     oldspin, newspin,   current_x, current_y,   nn_i, nn_j; // Lattice indices
    int root_x, root_y,     coord_x, coord_y,   pos_update_x, pos_update_y; // Coordinates for wrapping criteria
    int old_x, old_y,       new_x, new_y;       // More coordinates for wrapping criteria 
    int sp, cluster_size; // Counters
    double randnum; // Number to decide if add to cluster
    bool wrapping_crit = false;
    vector<pair<int, int> > stacker;
    stacker.reserve(LEN*LEN);
    stack<int> cluster_x;
    stack<int> cluster_y;
    
    // Pick a random spin (x then y coord) and make it the origin (i.e. "root/seed/(0,0)")
    i = rand() % LEN;
    j = rand() % LEN;
    spinSite* root_spin = getSpinSite(i,j);
    root_spin->AddToCluster();
    root_spin->Triangulate(0,0);
    cluster_x.push(i);
    cluster_y.push(j);
    stacker.push_back(make_pair(i, j));
    cluster_size = 1;

    sp = 1;
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

        for (int k = 0; k < NN_MAX; k++)
        {   
            
            nn_i = neighbor_table[current_x][current_y][k][0];
            nn_j = neighbor_table[current_x][current_y][k][1];

            spinSite* neighbor_spin = getSpinSite(nn_i,nn_j);
            bool neighbor_checked = neighbor_spin->checkStatus();
            if (wrapping_crit == false && neighbor_checked == true) { // Check for wrapping once

                old_x = neighbor_spin->getX();
                old_y = neighbor_spin->getY();
                int pos_update_x = coord_x, pos_update_y = coord_y;
                if (k == 0) pos_update_x -= 1;
                else if (k == 1) pos_update_x += 1;
                else if (k == 2) pos_update_y += 1;
                else if (k == 3) pos_update_y -= 1;
                else if (k == 4) pos_update_y += 2;
                else if (k == 5) pos_update_y -= 2;

                neighbor_spin->Triangulate(pos_update_x, pos_update_y);
                new_x = neighbor_spin->getX();
                new_y = neighbor_spin->getY();

                if (old_x != new_x || old_y != new_y) {
                    wrapping_crit = true;
                    wrap_counter++;
                }
            }   else if (neighbor_checked == false)  {
                randnum = static_cast<double>(rand()) / RAND_MAX;
                if ((neighbor_spin->getSpin() == oldspin && k <= 3 && randnum <= *padd1) ||
                    (neighbor_spin->getSpin() == newspin && k >= 4 && randnum <= *padd2)) 
                {
                    
                    int pos_update_x = coord_x, pos_update_y = coord_y;
                    if (k == 0) pos_update_x -= 1;
                    else if (k == 1) pos_update_x += 1;
                    else if (k == 2) pos_update_y += 1;
                    else if (k == 3) pos_update_y -= 1;
                    else if (k == 4) pos_update_y += 2;
                    else if (k == 5) pos_update_y -= 2;

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

    // For simpler percolation data, we are just interested in the non-wrapping cluster size.
    if (wrapping_crit == false)
    {
        avg_nowrap_cluster_size += cluster_size;
        nowrap_counter++;
    }
    avg_cluster_size += cluster_size;

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

void Lattice::doSweep(double *Beta)
{
    // int Ene, E1, E2,    Mag, M1, M2, M1_abs, Mag_abs, Mag_sq;
    double padd1, padd2;
    padd1 = 1 - exp(-2 * *Beta * J);
    padd2 = 1 - exp(-2 * *Beta * J * kappa);
    
    // E1 = E2 = M1 = M2 = M1_abs = 0;
    for (int j = 0; j < STEPS; j++) 
    {
        doStep(&padd1, &padd2); // Consider making one 'sweep' as a number of steps, where we choose it as after each spin has had one opportunity on average to flip
    }
        /*
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
        */
}

void Lattice::doWolffAlgo(double *Beta, fftw_plan p)
{
    // Only do the MC steps. Burn In is completed during initialization.
    wrap_counter = 0, nowrap_counter = 0;
    avg_cluster_size = 0, avg_nowrap_cluster_size = 0;
    // int num_sweeps = (*Beta <= 0.46) ? SWEEPS : SWEEPS / 4;
    int num_sweeps = (*Beta < .25) ? SWEEPS*10 : (*Beta < .46) ? SWEEPS : SWEEPS/4;
    double padd1 = 1 - exp(-2 * *Beta * J);
    double padd2 = 1 - exp(-2 * *Beta * J * kappa);
    for (int i = 0; i < num_sweeps; i++)
    {
        for (int j = 0; j < STEPS; j++)
        {
            doStep(&padd1, &padd2);
        }
    }
    doFFT(p);
}

void Lattice::updateTotalEnergy()
{
    int left_x, right_x, up_x, down_x;   // nearest neighbours
    int left_y, right_y, up_y, down_y;
    int up2_x, down2_x; // next-nearest neighbours
    int up2_y, down2_y;

    energy = 0;
    for (int i = 0; i < LEN; i++)
    {
        for (int j = 0; j < LEN; j++)
        {
            left_x = neighbor_table[i][j][0][0]; // I feel like I don't need to do all this...
            left_y = neighbor_table[i][j][0][1]; // This is hardcoded and not good for generalizing to any dimension
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

            int root = mid_s->getSpin();
            int L = left_s->getSpin();
            int R = right_s->getSpin();
            int U = up_s->getSpin();
            int D = down_s->getSpin();
            int U2 = up2_s->getSpin();
            int D2 = down2_s->getSpin();


            energy += -J * (root * (L + R + U + D));
            energy += J * kappa * (root * (U2 + D2));
            // energy += (mid_s->getSpin())*(-1*J) * (((left_s->getSpin()) + (right_s->getSpin()) \
            //                                     + (up_s->getSpin()) + (down_s->getSpin())) \
            //                                     - (kappa * ((up2_s->getSpin()) + (down2_s->getSpin()))));
        }  
    }   
    energy = (double)(energy / 2.0);
}

void Lattice::updateTotalMag()
{
    int i, j;
    int mag_ph = 0; // placeholder
    for (i = 0; i < LEN; i++){
        for (j = 0; j < LEN; j++){
            spinSite* site = getSpinSite(i,j);
            mag_ph += (site->getSpin());
        }
    }
    mag = mag_ph;
}

void Lattice::doFFT(fftw_plan p)
{
    // Prepare the FFT (typical FFTW implementation)
    double *in, *out;
    in = (double*) fftw_alloc_real(LEN);
    out = (double*) fftw_alloc_real(LEN);
    // fftw_plan p;
    // p = fftw_plan_r2r_1d(LEN, in, out, FFTW_R2HC, FFTW_MEASURE);
    
    // Prepare input array (slices of lattice)
    for (int i = 0; i < LEN; i++)
    {
        double slice_mag = 0;
        for (int j = 0; j < LEN; j++)
        {
            spinSite* site = getSpinSite(j,i); // Correct order of ij axes to see modulation
            slice_mag += (site->getSpin());
        }
        slice_mag /= (double)LEN;
        in[i] = slice_mag;
    }

    // Do FFT
    fftw_execute_r2r(p, in, out);

    // Get array of frequencies and corresponding weights
    int halfLEN = LEN/2;
    double data[halfLEN];
    double freqs[halfLEN];
    for (int i = 0; i < halfLEN; i++)
    {
        freqs[i] = (double)i/LEN; // I think this makes the right frequency values (equiv. to np.linspace)
        if (i == 0 || i == halfLEN-1){
            data[i] = sqrt(out[i]*out[i]);
        } else {
            data[i] = sqrt(out[i]*out[i] + out[LEN-i]*out[LEN-i]); // Weird indexing because of halfcomplex array
        }
    }

    // Get frequency with max. weight (store both the freq. and its weight as separate observables)
    double max = 0;
    int argmax = 0;
    for (int j = 0; j < halfLEN; j++)
    {
        if (data[j] > max){
            max = data[j];
            argmax = j;
        }
    }
    dom_freq = freqs[argmax];
    dom_amplitude = max;
}

