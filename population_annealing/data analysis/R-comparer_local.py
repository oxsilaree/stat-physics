import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

L = 32
kappa = 0.6
kappastr = f"{kappa:.2f}"
quantity = 11


R_vals = np.array([50,100,500,1000,2500,5000,7500,10000,15000,20000])
R_vals = np.array([100,500,2500])
colormap = plt.cm.brg
num_colors = len(R_vals)
colors = [colormap(i / num_colors) for i in range(num_colors)]
markers = [".","x"]
dfs = []

# Populate dataframes
info_name = f"/Users/shanekeiser/Documents/ANNNI/populationannealing/data/t/parameter_info_{kappastr}_kappa_{L}_L_{R_vals[0]}_R.csv"
info = np.loadtxt(info_name, dtype = str, delimiter = ',', skiprows = 1)
culling_frac = info[-1]
for R in R_vals:
    fname= f"/Users/shanekeiser/Documents/ANNNI/populationannealing/data/t/emcx_data_{kappastr}_kappa_{L}_L_{R}_R.csv"
    df = pd.read_csv(fname)
    dfs.append(df)

header = list(dfs[0])[quantity]

fig, ax = plt.subplots()

for i in range(len(R_vals)):
    if quantity == 18:
        ax.plot(dfs[i]["Beta"], -1*dfs[i]["Free Energy"]/dfs[i]["Beta"], marker = markers[i%2], color = colors[i], linewidth = 0, label = f"R = {R_vals[i]}")
    else:
        ax.plot(dfs[i]["Beta"], dfs[i][header], marker = markers[i%2], color = colors[i], linewidth = 0, label = f"R = {R_vals[i]}")

ax.set_title(f"Comparison of {header} across population sizes" + "\n" + fr"$\kappa$ = {kappastr}, $L$ = {L}, Culling fraction = {culling_frac}")
ax.set_xlabel(r"$\beta$")
ax.set_ylabel(f"{header}")
plt.legend()
plt.tight_layout()
plt.show()





