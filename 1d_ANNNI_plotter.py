"""
    This file reads the data from several data files (output from an older 1D_ANNNI_Cluster variant).

    Now, this doesn't really serve much purpose, since 1D_ANNNI_Annealed.py replaces it. This is because
    I changed the code such that all data is stored in 1 larger file. But this was a good exercise in using
    glob to read files from a directory.
"""

import glob
import matplotlib.pyplot as plt
import numpy as np
print(glob.glob('/Users/shanekeiser/Documents/Spring 2024/Machta Spring/1d_data/*ANNNI*.csv'))

filenames = sorted(glob.glob('/Users/shanekeiser/Documents/Spring 2024/Machta Spring/1d_data/*.csv'))
print(len(filenames))
Cvals = []
Xvals = []
i = 0
datasets = 25
for f in filenames:
    data = np.loadtxt(fname=f, delimiter=',')
    Cvals.append(data[0])
    Xvals.append(data[1])
    i += 1
    if ((i % 25) == 0):
        print(f'Getting closer! {i}')

print(filenames[0])
print(filenames[100])
print(Cvals)
print(Xvals)
Cvals = np.array(Cvals)
Cvals /= 100 # Divide by L
Cvals = Cvals.reshape((datasets,datasets))
print(len(Cvals))

Xvals = np.array(Xvals)
Xvals /= 100 # Divide by L
Xvals = Xvals.reshape((datasets,datasets))

T_vals = np.linspace(0,2,datasets)
k_vals = np.linspace(0,2,datasets)

plt.pcolor(T_vals,k_vals, Cvals, cmap = 'gray')
plt.title('Specific heat capacity per spin')
plt.ylabel('$\kappa$')
plt.xlabel('T')
plt.colorbar()
plt.show()


plt.pcolor(T_vals,k_vals, Xvals, cmap = 'gray')
plt.title('Mag. Susceptibility per spin')
plt.ylabel('$\kappa$')
plt.xlabel('T')
plt.colorbar()
plt.show()


