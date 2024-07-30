import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

# data = np.loadtxt("./data/emcx_data_0_kappa.csv", delimiter = ',', dtype = str)
data = np.loadtxt("/Users/shanekeiser/Downloads/emcx_data_0_kappa-16.csv", delimiter = ',', dtype = str) # L = 16, pop_size = 200
data = np.loadtxt("/Users/shanekeiser/Downloads/emcx_data_0_kappa-17.csv", delimiter = ',', dtype = str) # L = 16, pop_size = 1000
data = np.loadtxt("/Users/shanekeiser/Downloads/emcx_data_0_kappa-18.csv", delimiter = ',', dtype = str) # L = 64, pop_size = 1000
data = np.loadtxt("/Users/shanekeiser/Downloads/emcx_data_0_kappa-27.csv", delimiter = ',', dtype = str) # L = 16, pop_size = 500
data = np.loadtxt("/Users/shanekeiser/Downloads/emcx_data_2_kappa-2.csv", delimiter = ',', dtype = str)
data = np.loadtxt("/Users/shanekeiser/Downloads/emcx_data_0.5_kappa.csv", delimiter = ',', dtype = str)
data = np.loadtxt("/Users/shanekeiser/Downloads/emcx_data_1_kappa.csv", delimiter = ',', dtype = str)
data = np.loadtxt("/Users/shanekeiser/Downloads/emcx_data_1.5_kappa.csv", delimiter = ',', dtype = str)
trimmed_data = data[:,1:-1].astype(float)
print(data[:,0])
# print(trimmed_data[0])
B_vals = trimmed_data[0]
E_vals = trimmed_data[1]
E2_vals = trimmed_data[2]
M_vals = trimmed_data[3]
M2_vals = trimmed_data[4]
MA_vals = trimmed_data[5]
C_vals = trimmed_data[6]
X_vals = trimmed_data[7]
CS_vals = trimmed_data[8]
NWCS_vals = trimmed_data[9]
W_vals = trimmed_data[10]
print(B_vals)

def plotter(arr, header, opacity = 1, marker = '.', x = B_vals, linewidth = 0):
    data = arr
    temps = x
    plt.plot(temps, data, marker = marker, label = header, alpha = opacity, linewidth = linewidth)
    plt.xlabel(r"$\beta$")
    plt.axvline(1/2.2691853, color = 'black', alpha = 0.4)
    plt.xticks([0.3,1/2.2691853,0.6,0.9,1.2])
    plt.legend()

L=16
N=L**2
CULLING_FRAC=0.02
INIT_POP_SIZE=500
kappa = 1.5
titlestr = f"kappa = {kappa}, L = {L}, INIT_POP_SIZE = {INIT_POP_SIZE}, CULLING_FRAC = {CULLING_FRAC}"

### ENERGY
plotter(E_vals, "Energy")
plt.title(f"Energy per spin\n{titlestr}")
plt.show()

### ENERGY SQUARED
plotter(E2_vals, "Energy squared")
plt.title(f"Energy squared\n{titlestr}")
plt.show()

### SPECIFIC HEAT per spin
C = (E2_vals - 64*np.square(E_vals))*B_vals*B_vals
print(len(C_vals))
# C = C.astype(str)
# plotter(C, "Calculated in Python spec. heat", opacity = 0.5, marker = '*')
plotter(C_vals, "Calculated in C++", opacity = 0.5)
plt.show()

### MAGNETIZATION
plotter(M_vals, "Magnetization")
plotter(MA_vals, "Absolute magnetization")
plt.show()

### SQUARE MAGNETIZATION and CLUSTER SIZE
plotter(M2_vals, "Magnetization squared", marker = '.', linewidth = 1)
plotter(CS_vals, "Cluster size", marker = 'o', opacity = 0.5)
plt.title(f"Magnetization Squared vs. Cluster size\n{titlestr}")
plt.show()

### SUSCEPTIBILITY and NON-WRAPPING CLUSTER SIZE
X1 = (M2_vals - N*np.square(MA_vals))*B_vals
# plotter(X1, "Susceptibility (calculated in Python)", linewidth = 0.4)
plotter(X_vals, "Susceptibility (calculated in C++)", opacity = 0.5)
plotter(NWCS_vals * B_vals, r"Non-wrapping cluster size $\times \beta$")
# plt.plot(B_vals + 0.05, NWCS_vals * B_vals)
plt.title(f"Susceptibility plots \n{titlestr}")
plt.show()

### WRAPPING PERCENTAGE
plotter(E2_vals/(4*N), "Square energies / 4N")
plotter(M2_vals/N, "Square magnetizations / N")
# plotter(MA_vals, "Absolute magnetizations")
plotter(W_vals, "Fraction of clusters that wrap")
# plotter(CS_vals/N, "Cluster sizes / N")
plt.title(f'Wrapping probabilities comparison with other quantities\n{titlestr}')
plt.show()

### CLUSTER QUANTITIES

plotter(CS_vals, "Cluster sizes")
plotter(NWCS_vals, "Non-wrapping cluster sizes")
# plotter(W_vals, "Wrapping percentage")
plt.title(f"Comparison of wrapping vs. non-wrapping cluster sizes\n{titlestr}")
parent_axes = plt.gca()
axins = zoomed_inset_axes(parent_axes, 6, loc=5)
axins.plot(B_vals, CS_vals, '.')
axins.plot(B_vals, NWCS_vals, '.')
axins.set_xlim(0.28,0.41)
axins.set_ylim(0.02*N, 0.095*N)
for axis in ['top','bottom','left','right']:
    axins.spines[axis].set_linewidth(1)
    axins.spines[axis].set_color('k')
mark_inset(parent_axes, axins, loc1=2, loc2=4, fc="none", lw=1, ec='k')
# plt.draw()
plt.show()



fig, ax = plt.subplots(3,2)

fig.suptitle(f'Quantities (per spin)\n{titlestr}')
ax[0,0].plot(B_vals, E_vals, 'r.')
ax[0,0].set_xlabel(r'$\beta$')
ax[0,0].set_ylabel('Energy')
ax[0,0].set_xlim(left = np.min(B_vals))

ax[1,0].plot(B_vals, np.abs(MA_vals), 'r')
ax[1,0].set_xlabel(r'$\beta$')
ax[1,0].set_ylabel('(Absolute)\naverage Magnetization')
ax[1,0].set_xlim(left = np.min(B_vals))

ax[0,1].plot(B_vals, C_vals, 'r')
ax[0,1].set_xlabel(r'$\beta$')
ax[0,1].set_ylabel('Specific heat')
ax[0,1].set_xlim(left = np.min(B_vals))
ax[0,1].axvline(1/2.269, color = 'black')
# ax[0,1].text(2.3, 0.0, "T = 2.269", rotation = 270)

ax[1,1].plot(B_vals, X_vals/B_vals, 'r')
ax[1,1].set_xlabel(r'$\beta$')
ax[1,1].set_ylabel('Susceptibility')
ax[1,1].set_xlim(left = np.min(B_vals))
ax[1,1].axvline(1/2.269, color = 'black')
# ax[1,1].text(2.3, 0.0, "T = 2.269", rotation = 270)

ax[2,0].plot(B_vals, CS_vals, 'r')
ax[2,0].set_xlabel(r'$\beta$')
ax[2,0].set_ylabel('Average\ncluster size')
ax[2,0].set_xlim(left = np.min(B_vals))
ax[2,0].axvline(1/2.269, color = 'black')
# ax[2,0].text(2.3, 0.0, "T = 2.269", rotation = 270)

ax[2,1].plot(B_vals, NWCS_vals, 'r')
ax[2,1].set_xlabel(r'$\beta$')
ax[2,1].set_ylabel('Average non-wrapping\ncluster size')
ax[2,1].set_xlim(left = np.min(B_vals))
ax[2,1].axvline(1/2.269, color = 'black')
# ax[2,1].text(2.3, 0.0, "T = 2.269", rotation = 270)

plt.tight_layout()
plt.show()

print(len(B_vals))


multiplier = np.max(X_vals)/np.max(NWCS_vals)
plt.plot(B_vals,  X_vals, 'r', label = 'Susceptibility')
# plt.plot(B_vals, NWCS_vals, 'k', label = f'Non-wrap Cluster Size')
plt.plot(B_vals, NWCS_vals*B_vals, 'g', label = r'Non-wrap Cluster Size * $\beta$')
plt.plot(B_vals, NWCS_vals*multiplier, color = 'b', label = f'Non-wrap Cluster Size * max(susc)/max(NWCS)\n = {multiplier:.4f}')
plt.axvline(1/2.269, color = 'black', linestyle = 'dashed', alpha = 0.4)
plt.legend(loc = 'upper right')
plt.title("Susc. and Non-Wrapping Cluster Size Comparison\nL=32, init_pop_size = 250, CULLING_FRAC = 0.1")
plt.xlabel(r'$\beta$')
plt.tight_layout()
plt.xticks(np.array([0.3,1/2.2691853,0.6,0.9]))
plt.show()

N = 64
plt.plot(B_vals,M_vals, 'r', label = 'magnetization per spin')
plt.plot(B_vals,(CS_vals/N), 'b.', label = 'cluster size (normalized to M)')
plt.plot(B_vals,(CS_vals/N)**0.6, 'b.', label = 'cluster size (normalized to M, ^0.6)', alpha = 0.3)
plt.title('Magnetization and cluster size comparison\nL=32, init_pop_size = 250, CULLING_FRAC = 0.1')
plt.xlabel(r"$\beta$")
plt.legend(loc = 'lower right')
plt.show()

M2_vals = (X_vals/B_vals) + M_vals**2

plt.plot(B_vals,M2_vals*B_vals, 'r', label = 'magnetization per spin')
plt.plot(B_vals, X_vals, 'g', label = 'Susceptibility')
plt.plot(B_vals,(NWCS_vals*0.7), 'b.', label = 'non-wrapping cluster size (normalized to M)')
# plt.plot(B_vals,(CS_vals/N)**0.5, 'b.', label = 'cluster size (normalized to M, ^0.6)', alpha = 0.2)
plt.title('Magnetization and cluster size comparison\nL=32, init_pop_size = 250, CULLING_FRAC = 0.1')
plt.xlabel(r"$\beta$")
plt.legend(loc = 'lower right')
plt.show()

def multiplotter(fname = "nil"):
    if fname == "nil":
        print("No filename included. Exiting...")
        return (1)
    df = pd.read_csv(fname)
    df.head()



