import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

mode = "p"
kappa = 0
kappastr = f'{kappa:.2f}'
size = 16

method = "Undeclared"

if mode == "t":
    method = "Two-Replica with PA"
elif mode == "p":
    method = "Wolff with PA"

# fname = f"/Users/shanekeiser/Downloads/data/25-10-24/" + mode + f"/emcx_data_{kappastr}_kappa_{size}_L.csv"
# param_info = f"/Users/shanekeiser/Downloads/data/25-10-24/" + mode + f"/parameter_info_{kappastr}_kappa_{size}_L.csv"
fname = f"/Users/shanekeiser/Documents/ANNNI/populationannealing/data/" + mode + f"/emcx_data_{kappastr}_kappa_{size}_L.csv"
param_info = f"/Users/shanekeiser/Documents/ANNNI/populationannealing/data/" + mode + f"/parameter_info_{kappastr}_kappa_{size}_L.csv"


info = np.loadtxt(param_info, skiprows = 1, delimiter = ',', dtype = str)
df = pd.read_csv(fname)

pop_size = info[2]
# print(df.head)


beta_vals = df["Beta"]

# Cluster sizes
sizes = df[["Non-Wrapping Cluster Size", "X-Wrapping Cluster Size", "Z-Wrapping Cluster Size", "Both-Wrapping Cluster Size"]]
sizes.rename(columns={"Non-Wrapping Cluster Size": "xbar_zbar", "X-Wrapping Cluster Size": "x_zbar",
               "Z-Wrapping Cluster Size": "xbar_z", "Both-Wrapping Cluster Size": "x_z"}, inplace = True)

# Wrapping probabilities
probs = df[['Wrapping Prob. (Neither)','Wrapping Prob. (X-dir.)',
            'Wrapping Prob. (Z-dir.)','Wrapping Prob. (both dir.)']]
probs.rename(columns={'Wrapping Prob. (Neither)': "xbar_zbar",'Wrapping Prob. (X-dir.)': "x_zbar",
       'Wrapping Prob. (Z-dir.)': "xbar_z", 'Wrapping Prob. (both dir.)': "x_z"}, inplace = True)

wrap_options = {0 : "xbar_zbar",
                1 : "x_zbar",
                2 : "xbar_z",
                3 : "x_z",
                4 : "x_wrap",
                5 : "z_wrap"}


z_wrap_prob = probs["xbar_z"] + probs["x_z"]
x_wrap_prob = probs["x_zbar"] + probs["x_z"]
wrap_prob = probs["x_zbar"] + probs["x_z"] + probs["xbar_z"]
total_prob = probs["x_zbar"] + probs["x_z"] + probs["xbar_z"] + probs["xbar_zbar"]
plt.plot(beta_vals, z_wrap_prob, '.', label = r"All with Z wrapping ($xz + \bar{x}z$)")
plt.plot(beta_vals, x_wrap_prob, '.', label = r"All with X wrapping ($xz + x\bar{z}$)")
plt.plot(beta_vals, wrap_prob, '.', label = r"Any wrapping ($xz + \bar{x}z + x\bar{z}$)")
plt.plot(beta_vals, total_prob, '.', label = f"All possibilities\n" + r"($xz + \bar{x}z + x\bar{z} + \bar{x}\bar{z}$)")
plt.plot(beta_vals, probs["xbar_zbar"], '.', label = r"No wrapping ($\bar{x} \bar{z}$)")
plt.legend()
plt.title(f"Wrapping plots, {method}\nKappa = {kappastr}, lattice size = {size}, population size = {pop_size}")
plt.show()


plt.plot(beta_vals, probs["xbar_zbar"], '.', label = r"$\bar{x} \bar{z}$")
plt.plot(beta_vals, probs["xbar_z"], '.', label = r"$\bar{x} z$")
plt.plot(beta_vals, probs["x_zbar"], '.', label = r"$x\bar{z}$")
plt.plot(beta_vals, probs["x_z"], '.', label = r"$xz$")
plt.plot(beta_vals, total_prob, '.', label = f"All possibilities\n" + r"($xz + \bar{x}z + x\bar{z} + \bar{x}\bar{z}$)")
plt.legend()
plt.title(f"Disjoint wrapping plots, {method}\nKappa = {kappastr}, lattice size = {size}, population size = {pop_size}")
plt.show()
probs = probs.fillna(0)
sizes = sizes.fillna(0) ### NEED THIS?

z_wrap_sizes = probs["xbar_z"]*sizes["xbar_z"] + probs["x_z"]*sizes["x_z"]

# z_wrap_sizes = z_wrap_prob*(sizes["xbar_z"] + sizes["x_z"])

x_wrap_sizes = probs["x_zbar"]*sizes["x_zbar"] + probs["x_z"]*sizes["x_z"]

# x_wrap_sizes = x_wrap_prob*(sizes["x_zbar"] + sizes["x_z"])
wrap_sizes = probs["x_zbar"]*sizes["x_zbar"] + probs["x_z"]*sizes["x_z"] + probs["xbar_z"]*sizes["xbar_z"] 

total_sizes = probs["x_zbar"]*sizes["x_zbar"] + probs["x_z"]*sizes["x_z"] + probs["xbar_z"]*sizes["xbar_z"]  + probs["xbar_zbar"]*sizes["xbar_zbar"]


plt.plot(beta_vals, z_wrap_sizes, '.', label = "all clusters wrapping in Z")
plt.plot(beta_vals, x_wrap_sizes, '.', label = "all clusters wrapping in X")
plt.plot(beta_vals, wrap_sizes, '.', label = "all clusters wrapping")
plt.plot(beta_vals, total_sizes, '.', label = "all clusters")
plt.legend()
plt.xlim(min(beta_vals), max(beta_vals))
plt.title("Average wrapping sizes")
plt.show()

print(z_wrap_sizes)
plt.plot(beta_vals, sizes["xbar_z"], '.', label = "wrapping in only Z")
plt.plot(beta_vals, sizes["x_zbar"], '.', label = "wrapping in only X")
plt.plot(beta_vals, sizes["x_z"], '.', label = "wrapping in both dir")
plt.plot(beta_vals, sizes["xbar_zbar"], '.', label = "wrapping in neither dir")
plt.legend()
plt.title("Average wrapping sizes")
plt.show()

print(sizes["xbar_z"])