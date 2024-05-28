import numpy as np
import matplotlib.pyplot as plt

# Filenames subject to change
f1 = "/Users/shanekeiser/Documents/Spring 2024/Machta Spring/2D_ANNNI/2d_lattices/2d_ANNNI_annealing_visual_64kappahalf.csv"
f2 = "/Users/shanekeiser/Documents/Spring 2024/Machta Spring/2D_ANNNI/2d_lattices/2d_ANNNI_annealing_visual_64kappazero.csv"
f3 = "/Users/shanekeiser/Documents/Spring 2024/Machta Spring/2D_ANNNI/2d_lattices/2d_ANNNI_annealing_visual_64kappatwo.csv"
f4 = "/Users/shanekeiser/Documents/Spring 2024/Machta Spring/2D_ANNNI/2d_lattices/2d_ANNNI_annealing_visual.csv"

L = 64 # side length of lattice, to be updated depending on the size used in accompanying C++ script

def visualize(f,L):
    data = np.loadtxt(f, delimiter = ',', dtype = str)
    data = data[:-1]
    data = data.astype(int)
    data = data.reshape(9,L,L)
    return data

data1 = visualize(f1)

plt.imshow(data1[0,:,:])
plt.title("kappa = 1/2, initial")
plt.show()

plt.imshow(data1[-1,:,:])
plt.title("kappa = 1/2, end state")
plt.show()

data2 = visualize(f2)

plt.imshow(data2[0,:,:])
plt.title("kappa = 0, initial")
plt.show()

plt.imshow(data2[-1,:,:])
plt.title("kappa = 0, end state")
plt.show()

data3 = visualize(f3)

plt.imshow(data3[0,:,:])
plt.title("kappa = 2, end state")
plt.show()

plt.imshow(data3[-1,:,:])
plt.title("kappa = 2, end state")
plt.show()

### Subplots

fig, ax = plt.subplots(2,3)


ax[0,0].imshow(data1[0,:,:])
ax[0,0].set_title("$\kappa$ = 1/2, initial")
ax[0,0].axis('off')

ax[0,1].imshow(data2[0,:,:])
ax[0,1].set_title("$\kappa$ = 0, initial")
ax[0,1].axis('off')

ax[0,2].imshow(data3[0,:,:])
ax[0,2].set_title("$\kappa$ = 2, initial")
ax[0,2].axis('off')


ax[1,0].imshow(data1[-1,:,:])
ax[1,0].set_title("$\kappa$ = 1/2, end state")
ax[1,0].axis('off')

ax[1,1].imshow(data2[-1,:,:])
ax[1,1].set_title("$\kappa$ = 0, end state")
ax[1,1].axis('off')

ax[1,2].imshow(data3[-1,:,:])
ax[1,2].set_title("$\kappa$ = 2, end state")
ax[1,2].axis('off')

plt.tight_layout()
plt.show()


data4 = visualize(f4)


plt.imshow(data4[0,:,:])
# plt.title("kappa = 2, end state")
plt.show()

plt.imshow(data4[-1,:,:])
# plt.title("kappa = 2, end state")
plt.show()

