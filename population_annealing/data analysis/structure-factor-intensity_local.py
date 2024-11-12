import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

kappa = 0.6
kappastr = f"{kappa:.2f}"
size = 48
fname = f"/Users/shanekeiser/Documents/ANNNI/populationannealing/data/t/sfi_data_{kappastr}_kappa_{size}_L.csv"
fname_t= f"/Users/shanekeiser/Documents/ANNNI/populationannealing/data/t/emcx_data_{kappastr}_kappa_{size}_L.csv"

data = np.loadtxt(fname, delimiter = ',', dtype = str)[:,:-1]

df = pd.read_csv(fname_t)
beta = df["Beta"]
freqs = data[0,:]
freq_names = []

freq_nums = np.arange(0,int(size/4)+1,1)
for i in range(len(freqs)):
    freq_string = f"{freq_nums[i]}" + f"/{size}"
    freq_names.append(freq_string)
print(freq_names)
print(freq_nums)

data = np.array(data[1:,:], dtype = float)

markers = [".","x","^","+","v","o","*","1","2","3","4","8","s","p"]

for i in range(len(freqs)):
    plt.plot(beta, data[:-1,i], marker = markers[i], linewidth = 0, label = f"{freq_names[i]}")
plt.legend()
plt.title(f"Structure Factor Intensity plot for\nkappa = {kappa}, L = {size}")
plt.xlabel(r"$\beta$")
plt.ylabel("Counts")
plt.show()