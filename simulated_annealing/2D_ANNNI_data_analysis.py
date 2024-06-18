"""

  This script runs simple data analysis for the 2D ANNNI model. Particularly it gets
  Energy, Magnetization, Specific Heat and Susceptibility. I also attempt to fit interpolated 
  data to a general cosine function, to get the amplitude/wave number, as an order parameter.

"""


import glob
import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import rfft,rfftfreq
from scipy.interpolate import make_interp_spline, BSpline
from scipy.optimize import curve_fit


filenames = sorted(glob.glob('/Users/shanekeiser/Documents/Summer 2024/Research/Annealing/multirun/data_folders/data/*data*.csv'))
print(filenames)

L = 32
N = L**2
datasets = 100
n_kappas = datasets # Each dataset corresponds to a different kappa value
n_temps = 100

def Analyze(f,blocks): ### FOR E,M,C,X
    data = np.loadtxt(fname = f, delimiter = ',', dtype = str)
    data_avgs = np.empty(shape = (n_temps,5))
    for i in range(5):
        data_avgs[:,i] = Shortcut(data,i,blocks)  ### 8 blocks
    return data_avgs    ### returns [E,M,abs(M),C,X]

def Shortcut(arr,ind,blocks):
    if ind == None:
        vals = arr[1:-1].astype(float)
    else:
        vals = arr[ind,1:-1].astype(float)
    vals_avgs = np.zeros(int(len(vals)/blocks))
    for i in range(len(vals_avgs)):
        vals_avgs[i] = np.mean(vals[blocks*i:blocks*(i+1)])
    return vals_avgs

### Note that the files go down in kappa and down in temperature
### so we need to flip the order of both
data_vals = np.empty(shape = (5,datasets,n_temps))

b1 = 100 ### Number of blocks
for i in range(5):
    for j in range(datasets):
        data_vals[i,j,:] = Analyze(filenames[j],b1)[:,i]

# transpose makes temperature and kappa in the right order
E_vals = np.transpose(data_vals[0,:,:])
# E_vals[E_vals > 0] = None
# E_vals[E_vals < -640] = None
M_vals = np.transpose(data_vals[2,:,:])
C_vals = np.transpose(data_vals[3,:,:])
X_vals = np.transpose(data_vals[4,:,:])

print(1)
fig, ax = plt.subplots(2,2)
fig.suptitle(f"E,M,C,X plots for L = {L} (N = {N})")

im_1 = ax[0,0].imshow(E_vals, extent = [0.02,2,0.03,3], cmap = 'coolwarm_r', aspect = 'auto')
plt.colorbar(im_1,ax = ax[0,0])
ax[0,0].set_title("Energies")
ax[0,0].set_ylabel("T")
ax[0,0].set_xlabel("$\kappa$")

im_2 = ax[0,1].imshow(np.sqrt(M_vals), extent = [0.02,2,0.03,3], cmap = 'coolwarm_r', aspect = 'auto')
plt.colorbar(im_2,ax = ax[0,1])
ax[0,1].set_title("Absolute Magnetizations (square rooted)")
ax[0,1].set_ylabel("T")
ax[0,1].set_xlabel("$\kappa$")

im_3 = ax[1,0].imshow(C_vals/N, extent = [0.02,2,0.03,3], cmap = 'coolwarm_r', aspect = 'auto')
plt.colorbar(im_3,ax = ax[1,0])
ax[1,0].set_title("Spec. heat per spin")
ax[1,0].set_ylabel("T")
ax[1,0].set_xlabel("$\kappa$")

im_4 = ax[1,1].imshow(np.log(X_vals), extent = [0.02,2,0.03,3], cmap = 'coolwarm_r', aspect = 'auto')
plt.colorbar(im_4,ax = ax[1,1])
ax[1,1].set_title("Log Susceptibilities")
ax[1,1].set_ylabel("T")
ax[1,1].set_xlabel("$\kappa$")

plt.tight_layout()
plt.show()

# odd because there seems to be no variation in such things based on changing kappa?

# SPECIFIC HEAT VS. TEMPERATURE PLOTS
print(filenames[12])
print(filenames[3])
print(filenames[18])
print(filenames[24])
print(filenames[0])
kappa_half_data = np.flipud(Analyze(filenames[12], b1)) ### These need to be flipped since the temperature data goes as decreasing temperature
kappa_underhalf_data = np.flipud(Analyze(filenames[4],b1))
kappa_overhalf_data = np.flipud(Analyze(filenames[18],b1))
kappa_one_data = np.flipud(Analyze(filenames[24],b1))
kappa_zero_data = np.flipud(Analyze(filenames[0],b1))
T_vals = np.linspace(0.03,3,100)
plt.plot(T_vals, kappa_half_data[:,3],'.', label = "kappa = 0.52")
plt.plot(T_vals, kappa_underhalf_data[:,3],'.', label = "kappa = 0.20")
plt.plot(T_vals, kappa_overhalf_data[:,3],'.', label = "kappa = 0.76")
plt.legend()
plt.show()

#Interpolated plots (look nicer)
spl1 = make_interp_spline(T_vals, kappa_half_data[:,3], k=3)
spl2 = make_interp_spline(T_vals, kappa_overhalf_data[:,3], k=3)
spl3 = make_interp_spline(T_vals, kappa_underhalf_data[:,3], k=3)
spl4 = make_interp_spline(T_vals, kappa_one_data[:,3], k=3)
spl5 = make_interp_spline(T_vals, kappa_zero_data[:,3], k=3)
T_new = np.linspace(0.03,3,500)
smooth1 = spl1(T_new)
smooth2 = spl2(T_new)
smooth3 = spl3(T_new)
smooth4 = spl4(T_new)
smooth5 = spl5(T_new)

plt.plot(T_new,smooth5/N,color = "#FEDF00", label = r"$\kappa$ = 0.00")
plt.plot(T_new,smooth3/N,color = "#D3CD00", label = r"$\kappa$ = 0.20")
plt.plot(T_new,smooth1/N,color = "#00D333", label = r"$\kappa$ = 0.52",)
plt.plot(T_new,smooth2/N,color = "#00D3AD", label = r"$\kappa$ = 0.76")
plt.plot(T_new,smooth4/N,color = "#007DD3", label = r"$\kappa$ = 1.00")
plt.legend()
plt.title(f"Spec. heat per spin (L = {L})")
plt.xlabel("T")
plt.ylabel(r"spec. heat per spin")
plt.vlines(x = 2.269,ymin = 0,ymax = np.max(smooth3)/N,linestyle = 'solid', color = 'k')
plt.annotate(text = "T = 2.269", xy = [2.28,1.25], xytext = [3.2,1.2], arrowprops \
             = dict(arrowstyle='<-', color='black', linewidth=1))
plt.tight_layout()
plt.show()

# ORDER PARAMETER STUFF

def Modulator(z,A,k,phi): # Function for curve_fit
    return A*np.cos(2*np.pi*k*z + phi)

def AmplitudeIdentifier(input_arr,L): # This function analyzes one temperature/kappa data set
    nblocks = int(len(input_arr)/L)
    input_arr = input_arr.reshape(nblocks, L)
    output = np.empty(shape = (nblocks,3))
    for i in range(nblocks):
        arr = input_arr[i,:]

        z = np.arange(1,L+1,1)                 ##### Interpolates data to make an easier curve to fit cosine to
        z_new = np.linspace(z.min(), z.max(), 300)
        spl = make_interp_spline(z, arr, k = 3)
        y_smooth = spl(z_new)

        amps = rfft(arr)                        ##### Gets best frequency for cosine fitting guess
        freqs = rfftfreq(L, 1)
        max_freq = freqs[np.argmax(np.abs(amps))]

        pguesses = [0.5, max_freq, 1]           ###### Fits to cosine
        popt, pcov = curve_fit(Modulator, z_new, y_smooth, method = 'lm',maxfev = 16000, p0 = pguesses)
        output[i] = popt
    output[:,0] = np.abs(output[:,0])
    return output

nblocks = 100
L_1 = 32 ### Keep this separate from other L for now
mz_filenames = sorted(glob.glob('/Users/shanekeiser/Documents/Summer 2024/Research/Annealing/multirun/data_folders/data/*mz_vals*.csv'))

mz_files = np.empty(shape = (100,100*100*32))
for i in range(50):
    mz_files[i,:] = np.loadtxt(mz_filenames[i], delimiter = ',', dtype = str)[1:-1].astype(float)
mz_files = mz_files.reshape(100,100,100*32)


cosine_params = np.empty(shape = (100,100,100,3)) ### 50 kappa values, 100, T values, 100 blocks, 3 parameters

for i in range(100):
    for j in range(100):
        cosine_params[i,j,:,:] = AmplitudeIdentifier(mz_files[i,j,:], L = 32)
        if j%10 == 0:
            print(f"{i*100 + j} sets done!")

print(1)

cosine_param_averages = np.empty(shape = (100,100,3))
for i in range(100):
    for j in range(100):
        for k in range(3):
            cosine_param_averages[i,j,k] = np.average(cosine_params[i,j,:,k])

print(1)

# Plot relevant order parameter stuff
amplitudes = cosine_param_averages[:,:,0].T
print(np.min(amplitudes))
print(amplitudes.shape)
# amplitudes[amplitudes > 5] = 0
plottable_amps = np.log(amplitudes[:,:50] + 1)
plottable_amps[plottable_amps >= 3] = 1
plt.imshow(plottable_amps, extent = [0.02,1,0.03,3], aspect = 'auto')
plt.colorbar()
plt.title("Amplitude of modulated structure 'waveform' (L = 32)")
plt.xlabel(f"$\kappa$")
plt.ylabel("T")
plt.show()
wavenumbers = cosine_param_averages[:,:,1].T
plt.imshow(wavenumbers[:,:50], extent = [0.02,1,0.03,3], aspect = 'auto')
plt.colorbar()
plt.title("Wavenumber of modulated structure 'waveform' (L = 32)")
plt.xlabel(f"$\kappa$")
plt.ylabel("T")
plt.show()

bins = np.arange(0,1+(1/32),1/32) # Try a binned one since the phase should not be continuous, rather have specified values and transitions

wavenumbers_binned = np.digitize(wavenumbers,bins)/L
plt.imshow(wavenumbers_binned[:,:50], extent = [0.02,1,0.03,3], aspect = 'auto')
plt.title("Wavenumber of modulated structure 'waveform' (L = 32),\nbinned into multiples of 1/32 (1/L)")
plt.xlabel(f"$\kappa$")
plt.ylabel("T")
plt.colorbar()
plt.show()
       
### Find average CLUSTER SIZES

clusterfiles = sorted(glob.glob('/Users/shanekeiser/Documents/Summer 2024/Research/Annealing/multirun/data_folders/data/*cluster*.csv'))
cluster_data = np.empty(shape = (100,100))

for i in range(100):
    cluster_data[i,:] = np.loadtxt(fname = clusterfiles[i], delimiter = ',', dtype = str)[:-1].astype(float)

cluster_data = cluster_data.T
plt.imshow(cluster_data[:,:50], extent = [0.02,1,0.03,3], cmap = 'coolwarm_r', aspect = 'auto')
plt.title(f'Average cluster size, for L = {L} (N = {N})')
plt.ylabel('T')
plt.xlabel(r'$\kappa$')
plt.colorbar()
plt.show()


