import numpy as np
import matplotlib.pyplot as plt

data1 = np.loadtxt("/Users/shanekeiser/Documents/Spring 2024/Machta Spring/1d_data/lattices/1d_ANNNI_k4.000_T0.010visual.csv", delimiter = ',', dtype = str)
data1 = data1[:-1]
print(len(data1))
data2 = np.empty(len(data1))
for i in range(len(data1)):
    data2[i] = data1[i].astype(int)
print(data2.dtype)
print(data2.shape)


data2 = data2.reshape((16,100))


plt.imshow(data2)
plt.show()

data3 = np.loadtxt("/Users/shanekeiser/Documents/Spring 2024/Machta Spring/1d_data/lattices/1d_ANNNI_k5.000_T0.001visual.csv", delimiter = ',', dtype = str)
# data3 = data3[:-1]
data3 = data3.astype(int)
data3 = data3.reshape((16,100))
plt.imshow(data3)
plt.show()

data4 = np.loadtxt("/Users/shanekeiser/Documents/Spring 2024/Machta Spring/1d_data/lattices/1d_ANNNI_k0.000_T0.001visual.csv", delimiter = ',', dtype = str)
# data3 = data3[:-1]
print(len(data4))
data4 = data4.astype(int)
data4 = data4.reshape((17,100))
plt.imshow(data4)
plt.show()

data5 = np.loadtxt("/Users/shanekeiser/Documents/Spring 2024/Machta Spring/1d_data/lattices/1d_ANNNI_k0.500_T0.001visual.csv", delimiter = ',', dtype = str)
# data3 = data3[:-1]
data5 = data5.astype(int)
data5 = data5.reshape((17,100))
plt.imshow(data5)
plt.show()

data6 = np.loadtxt("/Users/shanekeiser/Documents/Spring 2024/Machta Spring/1d_data/lattices/1d_ANNNI_k2.000_T0.001visual.csv", delimiter = ',', dtype = str)
# data3 = data3[:-1]
data6 = data6.astype(int)
data6 = data6.reshape((17,100))
plt.imshow(data6)
plt.show()


fig, ax = plt.subplots(3, 1, figsize = (9,6))

ax[0].imshow(data6,cmap = 'cividis')
ax[0].set_title(f"$\kappa$ = 2, $T$ = 0.001")
ax[0].set_ylabel("Epoch number")

ax[1].imshow(data5,cmap = 'cividis')
ax[1].set_title(f"$\kappa$ = 0.5, $T$ = 0.001")
ax[1].set_ylabel("Epoch number")

ax[2].imshow(data4,cmap = 'cividis')
ax[2].set_title(f"$\kappa$ = 0, $T$ = 0.001")
ax[2].set_ylabel("Epoch number")
# #5D1E05
for a in np.arange(0.5,15.5,1):
    ax[0].axhline(y=a, color = 'grey', linestyle = '-')
    ax[1].axhline(y=a, color = 'grey', linestyle = '-')
    ax[2].axhline(y=a, color = 'grey', linestyle = '-')

for b in np.arange(0.5,99.5,1):
    ax[0].axvline(x=b, color = 'grey', linestyle = '-')
    ax[1].axvline(x=b, color = 'grey', linestyle = '-')
    ax[2].axvline(x=b, color = 'grey', linestyle = '-')

plt.xlabel("Spin position")
plt.tight_layout()
plt.show()





f1 = "/Users/shanekeiser/Documents/Spring 2024/Machta Spring/1d_data/lattices/1d_ANNNI_k0.000_T3.000visual.csv"
f2 = "/Users/shanekeiser/Documents/Spring 2024/Machta Spring/1d_data/lattices/1d_ANNNI_k0.500_T3.000visual.csv"
f3 = "/Users/shanekeiser/Documents/Spring 2024/Machta Spring/1d_data/lattices/1d_ANNNI_k2.000_T3.000visual.csv"

def visualize(f):
    data = np.loadtxt(f, delimiter = ',', dtype = str)
    data = data.astype(int)
    data = data.reshape(17,100)
    return data

data7 = visualize(f1)
data8 = visualize(f2)
data9 = visualize(f3)


fig, ax = plt.subplots(3, 1, figsize = (9,6))

ax[0].imshow(data9,cmap = 'cividis')
ax[0].set_title(f"$\kappa$ = 2, $T$ = 3")
ax[0].set_ylabel("Epoch number")

ax[1].imshow(data8,cmap = 'cividis')
ax[1].set_title(f"$\kappa$ = 0.5, $T$ = 3")
ax[1].set_ylabel("Epoch number")

ax[2].imshow(data7,cmap = 'cividis')
ax[2].set_title(f"$\kappa$ = 0, $T$ = 3")
ax[2].set_ylabel("Epoch number")
# #5D1E05
for a in np.arange(0.5,15.5,1):
    ax[0].axhline(y=a, color = 'grey', linestyle = '-')
    ax[1].axhline(y=a, color = 'grey', linestyle = '-')
    ax[2].axhline(y=a, color = 'grey', linestyle = '-')

for b in np.arange(0.5,99.5,1):
    ax[0].axvline(x=b, color = 'grey', linestyle = '-')
    ax[1].axvline(x=b, color = 'grey', linestyle = '-')
    ax[2].axvline(x=b, color = 'grey', linestyle = '-')

plt.xlabel("Spin position")
plt.tight_layout()
plt.show()

f4 = "/Users/shanekeiser/Documents/Spring 2024/Machta Spring/1d_data/lattices/1d_ANNNI_k100.0_T0.001visual.csv"
data10 = visualize(f4)

plt.imshow(data10)
plt.show()


f5 = "/Users/shanekeiser/Documents/Spring 2024/Machta Spring/1d_data/lattices/1d_ANNNI_k10000_T0.001visual.csv"
data11 = visualize(f5)

plt.imshow(data11)
plt.show()


data12 = np.loadtxt("/Users/shanekeiser/Documents/Spring 2024/Machta Spring/1d_data/lattices/1d_ANNNI_k1.868_T0.250visual.csv", delimiter = ',', dtype = str)
# data3 = data3[:-1]
data12 = data12.astype(int)
data12 = data12.reshape((17,100))
plt.imshow(data12)
plt.show()
