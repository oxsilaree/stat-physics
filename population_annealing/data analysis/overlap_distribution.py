import numpy as np
import matplotlib.pyplot as plt
import pandas as pd



kappa = 1
kappastr = f'{kappa:.2f}'
size = 32

fname = f"/Users/shanekeiser/Documents/ANNNI/populationannealing/data/t/overlap_distribution_{kappastr}_kappa_{size}_L.csv"
param_info = f"/Users/shanekeiser/Documents/ANNNI/populationannealing/data/t/parameter_info_{kappastr}_kappa_{size}_L.csv"
emcx_data = f"/Users/shanekeiser/Documents/ANNNI/populationannealing/data/t/emcx_data_{kappastr}_kappa_{size}_L.csv"

df = pd.read_csv(emcx_data)

betas = df["Beta"]


# df = pd.read_csv(fname, sep = ',', header = None).T

data = np.loadtxt(fname, delimiter = ',', usecols = np.arange(0,2401,1))
print(data.shape)

temp_step = 255
print(betas[temp_step])
plt.hist(data[temp_step,:], bins = np.linspace(-1,1,200))
plt.xlabel("Overlap (normalized to lattice size)")
plt.ylabel("Counts")
plt.title(rf"Overlap distribution for $\beta$ = {betas[temp_step]}, $\kappa$ = {kappa}, L = {size}")
plt.show()