import numpy as np
import matplotlib.pyplot as plt
import pandas as pd



kappa = 0.55
kappastr = f'{kappa:.2f}'
size = 16

fname = f"/Users/shanekeiser/Documents/ANNNI/populationannealing/data/t/overlap_distribution_{kappastr}_kappa_{size}_L.csv"
param_info = f"/Users/shanekeiser/Documents/ANNNI/populationannealing/data/t/parameter_info_{kappastr}_kappa_{size}_L.csv"
emcx_data = f"/Users/shanekeiser/Documents/ANNNI/populationannealing/data/t/emcx_data_{kappastr}_kappa_{size}_L.csv"

df = pd.read_csv(emcx_data)

betas = df["Beta"]
info = np.loadtxt(param_info, delimiter = ',', dtype = str)

pop_size = info[1,2]

# df = pd.read_csv(fname, sep = ',', header = None).T

data = pd.read_csv(fname, delimiter = ',', names = range(5400))

temp_step = 100
print(betas[temp_step])
plt.hist(data[temp_step], bins = np.linspace(-1,1,100))
plt.xlabel("Overlap (normalized to lattice size)")
plt.ylabel("Counts")
plt.title(rf"Overlap distribution for $\beta$ = {betas[temp_step]}, $\kappa$ = {kappa}, L = {size}, R = {pop_size}")
plt.show()