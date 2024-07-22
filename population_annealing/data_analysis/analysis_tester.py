import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("./data/emcx_data_0_kappa.csv", delimiter = ',', dtype = str)

trimmed_data = data[:,1:-1].astype(float)
print(trimmed_data)

T_vals = 1/trimmed_data[0]
E_vals = trimmed_data[1]
M_vals = trimmed_data[2]
C_vals = trimmed_data[3]
X_vals = trimmed_data[4]

fig, ax = plt.subplots(2,2)

fig.suptitle('Quantities per spin')
ax[0,0].plot(T_vals, E_vals, 'r.')
ax[0,0].set_xlabel('T')
ax[0,0].set_ylabel('Energy')
ax[0,0].set_xlim(np.min(T_vals),4)

ax[1,0].plot(T_vals, np.abs(M_vals), 'r')
ax[1,0].set_xlabel('T')
ax[1,0].set_ylabel('(Absolute) average Magnetization')
ax[1,0].set_xlim(np.min(T_vals),4)

ax[0,1].plot(T_vals, C_vals, 'r')
ax[0,1].set_xlabel('T')
ax[0,1].set_ylabel('Specific heat')
ax[0,1].set_xlim(np.min(T_vals),4)
ax[0,1].axvline(2.269, color = 'black')
ax[0,1].text(2.3, 0.0, "T = 2.269", rotation = 270)

ax[1,1].plot(T_vals, X_vals, 'r')
ax[1,1].set_xlabel('T')
ax[1,1].set_ylabel('Susceptibility')
ax[1,1].set_xlim(np.min(T_vals),4)
ax[1,1].axvline(2.269, color = 'black')
ax[1,1].text(2.3, 0.0, "T = 2.269", rotation = 270)


plt.show()

print(len(T_vals))


