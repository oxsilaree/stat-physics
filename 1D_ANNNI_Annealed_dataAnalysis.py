"""
    This replaced 1d_ANNNI_plotter. This code accompanies 1D_ANNNI_Cluster_simAnnealing.cpp, which outputs 1 data 
    file.

    This file plots energy, magnetzn, spec. heat and suscepty after extracting that information from 1 data file.

    I leave some notes at the end for selecting a better 'order' parameter, but this is better saved for the 2D 
    ANNNI simulations.
"""


import glob
import matplotlib.pyplot as plt
import numpy as np
file = '/Users/shanekeiser/Documents/Spring 2024/Machta Spring/1d_data/1d_ANNNI_annealing.csv'

data = np.loadtxt(fname = file, delimiter = ',', dtype = str)
data = data[:,1:-1].astype(np.float64)

E = data[0,:]
M = data[1,:]
C = data[2,:]
X = data[3,:]

print(np.max(X))


runs = 50
# We have constant kappa for 10 runs, with T going down in increments of 0.2, and then kappa goes down by 0.2 each time too

def process(arr):
    arr_avgs = []
    for i in range(0,runs**2):
        j = 8*i
        arr_avgs.append(np.mean(arr[j:j+7]))
    arr_avgs = np.array(arr_avgs)
    arr_avgs = arr_avgs.reshape((runs,runs))

    arr_avgs = np.flip(arr_avgs,1)
    return arr_avgs

E_avgs = process(E)
M_avgs = process(M)
C_avgs = process(C)
X_avgs = process(X)
print(np.min(X_avgs))


X_avgs[X_avgs <= 1] = 1


print(np.min(E))
print(np.min(M))


fig, ax = plt.subplots(2,2, figsize = (8,8))

Eplot = ax[0,0].imshow(E_avgs, cmap= 'PuRd', interpolation = 'none', extent = [0,2,0,2])
ax[0,0].set_title("Energy values")
ax[0,0].set_ylabel("kappa values")
ax[0,0].set_xlabel("T values")
fig.colorbar(Eplot, ax=ax[0,0], shrink=0.7)


Mplot = ax[1,0].imshow(np.abs(M_avgs), cmap= 'PuRd', interpolation = 'none', extent = [0,2,0,2])
ax[1,0].set_title("Magnetization values")
ax[1,0].set_ylabel("kappa values")
ax[1,0].set_xlabel("T values")
fig.colorbar(Mplot, ax=ax[1,0], shrink=0.7)


Cplot = ax[0,1].imshow(C_avgs, cmap= 'PuRd', interpolation = 'none', extent = [0,2,0,2])
ax[0,1].set_title("Specific heat values")
ax[0,1].set_ylabel("kappa values")
ax[0,1].set_xlabel("T values")
fig.colorbar(Cplot, ax=ax[0,1], shrink=0.7)

Xplot = ax[1,1].imshow(np.log(X_avgs), cmap= 'PuRd', interpolation = 'none', extent = [0,2,0,2])
ax[1,1].set_title("Susceptibility values")
ax[1,1].set_ylabel("kappa values")
ax[1,1].set_xlabel("T values")
fig.colorbar(Xplot, ax=ax[1,1], shrink=0.7)

plt.tight_layout()
plt.show()


counter = 0
for i in range(len(X)):
    if X[i] == 0:
        counter += 1
print(counter)
print(len(X))

"""
Improved estimator:

in real life, such as in metropolis, the spins only flip one by one, so mag. is found 
for each phase


mag. associated with size of largest cluster

Only take absolute values of M for its data and plotting but still use regular values of M for X.

"""
