import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob


kappa = 0.6
L = 96
kappastr = f"{kappa:.2f}"
R = 1600
quantity = 24
mode = "t"
modestr = "Two-Replica_Method" if (mode == "t") else "Wolff_Method"

dirpath = f"/Users/shanekeiser/Downloads/production-run/{modestr}/{kappastr}_kappa/{L}_L/{R}_R"
filepattern = f"emcx_data_*"
files = glob.glob(f"{dirpath}/{filepattern}")

# info = np.loadtxt(info_name, delimiter = ',', dtype = str, skiprows=1)
r = np.random.randint(0, len(files)-1)
r1 = np.random.randint(0, len(files)-1)
df_t = pd.read_csv(files[r])
df_p = pd.read_csv(files[r1])
header = list(df_t)[quantity]


# text_box = f"kappa = {info[0]},\nL = {info[1]},\nINIT_POP_SIZE = {info[2]},\nCULLING_FRAC = {info[3]}"
fig, ax = plt.subplots()
ax.plot(df_t["Beta"], df_t[header], marker = '.', color = 'r', label = "Two replica", linewidth = 0)
ax.plot(df_p["Beta"], df_p[header], marker = '.', color = 'b', label = "Wolff", linewidth = 0)
ax.set_title(f"Comparison of {header}\nfor Two-Replica and Wolff algorithms")
ax.set_title(f"{header}\nfor Two-Replica")
ax.set_xlabel(r"Inverse temperature $\beta$")
ax.set_ylabel(f"{header}")
# ax.legend()
# ax.legend()
# ax.text(x = 0.5, y = 0.8, s = text_box,
#          fontsize=9, color="black", ha="left", va="center",
#          bbox=dict(facecolor="beige", edgecolor="black", boxstyle="round,pad=0.2"),
#          transform = ax.transAxes)
plt.show()

"""
betaF_t = -1*df_t["Free Energy"]/df_t["Beta"]
betaF_p = -1*df_p["Free Energy"]/df_p["Beta"]
print(betaF_t)
fig, ax = plt.subplots()
ax.plot(df_t["Beta"], betaF_t, marker = '.', color = 'r', label = "Two replica", linewidth = 0)
ax.plot(df_p["Beta"], betaF_p, marker = '.', color = 'b', label = "Wolff", linewidth = 0)
# ax.set_title(r"Comparison of free energy $F_k = - \frac{1}{\beta_k} \sum_{j=K}^{k+1} \ln Q(\beta_j, \beta_{j-1})$" + f"\n(for Two-Replica and Wolff algorithms)")
ax.set_title(r"Comparison of change in free energy $\Delta F_k = - \frac{1}{\beta_k} \ln Q(\beta_k, \beta_{k-1})$" + f"\n(for Two-Replica and Wolff algorithms)")
ax.legend()
ax.set_xlabel(r"Inverse temperature $\beta$")
ax.set_ylabel(r"$\Delta F$")
ax.text(x = 0.4, y = 0.5, s = text_box,
         fontsize=9, color="black", ha="left", va="center",
         bbox=dict(facecolor="beige", edgecolor="black", boxstyle="round,pad=0.2"),
         transform = ax.transAxes)
plt.tight_layout()
plt.show()
"""
