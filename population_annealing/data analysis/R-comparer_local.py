import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob

L = 48
kappa = 0.6
kappastr = f"{kappa:.2f}"
quantity = 6


R_vals = np.array([16000])
R = R_vals[0]
# R_vals = np.array([100,500,2500])
colormap = plt.cm.brg
num_colors = len(R_vals)
colors = [colormap(i / num_colors) for i in range(num_colors)]
markers = [".","x"]
dfs = []

# Populate dataframes
dirpath = f"/Users/shanekeiser/Documents/ANNNI/populationannealing/production-run/Two-Replica_Method/{kappastr}_kappa/{L}_L/{R}_R"
dirpath = f"/Users/shanekeiser/Downloads/production-run/Two-Replica_Method/0.60_kappa/48_L/{R}_R"
filepattern = f"emcx_data_*"
files = glob.glob(f"{dirpath}/{filepattern}")
info_name = f"parameter_info*"
print(files)
# info = np.loadtxt(dirpath + info_name, dtype = str, delimiter = ',', skiprows = 1)
# culling_frac = info[-1]
for i in range(len(R_vals)):
    df = pd.read_csv(files[i])
    dfs.append(df)

header = list(dfs[0])[quantity]

fig, ax = plt.subplots()

for i in range(len(R_vals)):
    if quantity == 18:
        ax.plot(dfs[i]["Beta"], -1*dfs[i]["Free Energy"]/dfs[i]["Beta"], marker = markers[i%2], color = colors[i], linewidth = 0, label = f"R = {R_vals[i]}")
    else:
        ax.plot(dfs[i]["Beta"], dfs[i][header], marker = markers[i%2], color = colors[i], linewidth = 0, label = f"R = {R_vals[i]}")

ax.set_title(f"Comparison of {header} across population sizes" + "\n" + fr"$\kappa$ = {kappastr}, $L$ = {L}")#, Culling fraction = {0.05}")
ax.set_title(f"{header} for R = 1600, S = 100" + "\n" + fr"$\kappa$ = {kappastr}, $L$ = {L}")#, Culling fraction = {0.05}")
ax.set_xlabel(r"$\beta$")
ax.set_ylabel(f"{header}")
# plt.legend()
ax.text(x = 0.15, y = 1.6, s = r"$\Delta \beta = 0.005$", bbox = dict(facecolor='yellow', alpha=0.5))
ax.text(x = 0.65, y = 1.6, s = r"$\Delta \beta = 0.0005$",bbox = dict(facecolor='yellow', alpha=0.5))
ax.vlines(x = [0.5], ymin = 0, ymax = 2, linestyle = 'dashed', alpha = 0.5)
plt.tight_layout()
plt.show()





