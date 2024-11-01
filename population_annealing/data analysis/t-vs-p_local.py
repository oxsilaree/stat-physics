import numpy as np
import matplotlib.pyplot as plt
import pandas as pd



kappa = 0.6
size = 16
kappastr = f"{kappa:.2f}"
quantity = 11

fname_t= f"/Users/shanekeiser/Documents/ANNNI/populationannealing/data/t/emcx_data_{kappastr}_kappa_{size}_L.csv"
fname_p= f"/Users/shanekeiser/Documents/ANNNI/populationannealing/data/p/emcx_data_{kappastr}_kappa_{size}_L.csv"
info_name = f"/Users/shanekeiser/Documents/ANNNI/populationannealing/data/t/parameter_info_{kappastr}_kappa_{size}_L.csv" # just use 't' for simplicity

info = np.loadtxt(info_name, delimiter = ',', dtype = str, skiprows=1)
df_t = pd.read_csv(fname_t)
df_p = pd.read_csv(fname_p)
header = list(df_t)[quantity]
text_box = f"kappa = {info[0]},\nL = {info[1]},\nINIT_POP_SIZE = {info[2]},\nCULLING_FRAC = {info[3]}"
fig, ax = plt.subplots()
ax.plot(df_t["Beta"], df_t[header], marker = '.', color = 'r', label = "Two replica", linewidth = 0)
ax.plot(df_p["Beta"], df_p[header], marker = '.', color = 'b', label = "Wolff", linewidth = 0)
ax.set_title(f"Comparison of {header}\nfor Two-Replica and Wolff algorithms")
ax.legend()
ax.text(x = 0.5, y = 0.8, s = text_box,
         fontsize=9, color="black", ha="left", va="center",
         bbox=dict(facecolor="beige", edgecolor="black", boxstyle="round,pad=0.2"),
         transform = ax.transAxes)
plt.show()

