import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob
from scipy.interpolate import make_interp_spline

N_steps = 1000
kappa = 0.6
kappastr = f"{kappa:.2f}"
L = 32
R = 1000 
mode = "t"
modestr = "Two-Replica_Method" if (mode == "t") else "Wolff_Method"

dirpath = f"/Users/shanekeiser/Downloads/production-run/{modestr}/{kappastr}_kappa/{L}_L/{R}_R"
filepattern = f"emcx_data_*"
files = glob.glob(f"{dirpath}/{filepattern}")


def getObs(fname, qty_number):
    df = pd.read_csv(fname)
    betas = df["Beta"]
    header = list(df)[qty_number]
    observable = df[header]
    obs_interp = make_interp_spline(betas,observable,k=3)

    beta_model = np.linspace(0.1,2,N_steps)
    obs_model = obs_interp(beta_model)
    return obs_model

def getlastFE(fname):
    df = pd.read_csv(fname)
    betas = df["Beta"]
    FE = -1*df["Free Energy"]/df["Beta"] # need to manipulate F a little, based on what is calculated in C++
    
    FE_interp = make_interp_spline(betas,FE,k=3)
    beta_model = np.linspace(0.1,2.5,N_steps)
    FE_model = FE_interp(beta_model)
    return FE_model[-200]

beta_model = np.linspace(0.1,2.5,N_steps)
obs_array = np.empty(shape = (len(files),N_steps))
FE_array = np.empty(len(files))
quantity = 11
for i in range(len(files)):
    obs_array[i] = getObs(files[i],quantity)
    FE_array[i] = getlastFE(files[i])

obs_means = np.empty(N_steps)
obs_stds = np.empty(N_steps)
for i in range(N_steps):
    obs_means[i] = np.mean(obs_array[:,i])
    obs_stds[i] = np.std(obs_array[:,i])

FE_mean = np.mean(FE_array)
FE_std = np.std(FE_array)
# print(f"<F> = {FE_mean} +/- {FE_std}")

# plt.errorbar(beta_model, obs_means, obs_stds)
# plt.show()


def getErrorBarPlot(kappa,L,R,mode, quantity):
    modestr = "Two-Replica_Method" if (mode == "t") else "Wolff_Method"
    dirpath = f"/Users/shanekeiser/Downloads/production-run/{modestr}/{kappastr}_kappa/{L}_L/{R}_R"
    filepattern = f"emcx_data_*"
    files = glob.glob(f"{dirpath}/{filepattern}")
    # beta_model = np.linspace(0.1,2.5,N_steps)
    obs_array = np.empty(shape = (len(files),N_steps))
    FE_array = np.empty(len(files))
    for i in range(len(files)):
        obs_array[i] = getObs(files[i],quantity)
        FE_array[i] = getlastFE(files[i])

    obs_means = np.empty(N_steps)
    obs_stds = np.empty(N_steps)
    for i in range(N_steps):
        obs_means[i] = np.mean(obs_array[:,i])
        obs_stds[i] = np.std(obs_array[:,i])

    return obs_means, obs_stds

def compareQuantity(kappa,L,R,quantity):
    kappastr = f"{kappa:.2f}"
    dirpath = f"/Users/shanekeiser/Downloads/production-run/Two-Replica_Method/{kappastr}_kappa/{L}_L/{R}_R"
    filepattern = f"emcx_data_*"
    filepattern_2 = f"param_info_*"
    files = glob.glob(f"{dirpath}/{filepattern}")
    df = pd.read_csv(files[0])
    # p_info = np.loadtxt(glob.glob(f"{dirpath}/{filepattern_2}")[0])
    header = list(df)[quantity]

    t_means, t_stds = getErrorBarPlot(kappa = kappa,L = L, R = R,mode = "t", quantity = quantity)
    p_means, p_stds = getErrorBarPlot(kappa = kappa,L = L, R = R,mode = "p", quantity = quantity)
    beta_vals = np.linspace(0.1,2.5,N_steps)
    
    
    plt.errorbar(beta_vals, t_means, t_stds, marker = '.', markersize = 0.5, color = 'b', label = "Two Replica", linewidth = 0)
    plt.errorbar(beta_vals, p_means, p_stds, marker = '.', markersize = 0.5, color = 'r', label = "Wolff", linewidth = 0)
    plt.errorbar(beta_vals, t_means, t_stds, color = 'b', linewidth = 0.5, alpha = 0.5)
    plt.errorbar(beta_vals, p_means, p_stds, color = 'r', linewidth = 0.5, alpha = 0.5)
    plt.legend()
    plt.title(f"Average over {len(files)} runs for {header}\n"\
              + f"L = {L}, R = {R}, " + r"$\kappa$" + f" = {kappa}")
    plt.show()


def getMeanFE(kappa,L,R,mode):
    kappastr = f"{kappa:.2f}"
    modestr = "Two-Replica_Method" if (mode == "t") else "Wolff_Method"
    dirpath = f"/Users/shanekeiser/Downloads/production-run/{modestr}/{kappastr}_kappa/{L}_L/{R}_R"
    filepattern = f"emcx_data_*"
    files = glob.glob(f"{dirpath}/{filepattern}")
    # beta_model = np.linspace(0.1,2.5,N_steps)
    FE_array = np.empty(len(files))
    for i in range(len(files)):
        FE_array[i] = getlastFE(files[i])/(L*L)


    FE_mean = np.mean(FE_array)
    FE_std = np.std(FE_array)
    # print(f"<F> = {FE_mean} +/- {FE_std}")

    return FE_mean, FE_std


""" FIX 
        THIS 
            FELLA"""
### Let's start from the beginning, shall we?

N_steps = 1000
beta_model = np.linspace(0.1,2,N_steps)


dirpath = f"/Users/shanekeiser/Downloads/production-run/{modestr}/{kappastr}_kappa/{L}_L/{R}_R"
filepattern = f"emcx_data_*"
files = glob.glob(f"{dirpath}/{filepattern}")

def getUnbiasedEstimate(kappa,L,R,mode,quantity):
    # This function returns the unbiased estimator for an observable, as well as its standard error which is found using the bootstrap method.
    kappastr = f"{kappa:.2f}"
    modestr = "Two-Replica_Method" if (mode == "t") else "Wolff_Method"
    dirpath = f"/Users/shanekeiser/Downloads/production-run/{modestr}/{kappastr}_kappa/{L}_L/{R}_R"
    filepattern = f"emcx_data_*"
    files = glob.glob(f"{dirpath}/{filepattern}")
    print(f"No. of datasets = {len(files)}")
    beta_model = np.linspace(0.1,2,N_steps)   
    df_example = pd.read_csv(files[0])
    header = list(df_example)[quantity]
    if L == 96:
        tail = 18
    else:
        tail = 4

    if quantity != 18 and quantity != 6:

        obs_matrix = np.zeros(shape = (len(files), N_steps))
        FE_weights = np.zeros(shape = (len(files), N_steps))
        FE_weight_denoms = np.zeros(shape = N_steps)
        FE_interp = np.zeros(shape = (len(files), N_steps))

        for i in range(len(files)):
            df = pd.read_csv(files[i])
            df.drop(df.tail(tail).index,inplace=True)
            betas = np.array(df["Beta"])
            observable = df[header]
            FE = df["Free Energy"] ### DIVIDE BY NUMBER OF SPINS TO MAKE NUMBERS MANAGEABLE --> Is this even OK to do? probably yes just to get weightings
            # print(FE)
            obs_spline = make_interp_spline(betas, observable, k=3)
            FE_spline = make_interp_spline(betas, FE, k=3)
            obs_interp = obs_spline(beta_model)
            obs_matrix[i,:] = obs_interp
            FE_interp[i,:] = FE_spline(beta_model)
            for j in range(N_steps):
                FE_weights[i,j] = np.exp(FE_interp[i,j] - max(FE_interp[:,j]))
                FE_weight_denoms[j] += FE_weights[i,j]
        # print(FE_weights[50,:])
    
        for i in range(N_steps):
            FE_weights[:,i] /= FE_weight_denoms[i]
        
        unbiased_obs = np.zeros(N_steps)
        for i in range(len(files)):
            unbiased_obs += obs_matrix[i,:] * FE_weights[i,:]
        
    if quantity == 6:

        ene_matrix = np.zeros(shape = (len(files), N_steps))
        ene_sq_matrix = np.zeros(shape = (len(files), N_steps))
        obs_matrix = np.zeros(shape = (len(files), N_steps))
        FE_weights = np.zeros(shape = (len(files), N_steps))
        FE_weight_denoms = np.zeros(shape = N_steps)
        FE_interp = np.zeros(shape = (len(files), N_steps))
        for i in range(len(files)):
            df = pd.read_csv(files[i])
            df.drop(df.tail(tail).index,inplace=True)
            betas = np.array(df["Beta"])
            energy = df["Energy"]
            energy_sq = df["Energy Squared"]
            FE = df["Free Energy"] ### DIVIDE BY NUMBER OF SPINS TO MAKE NUMBERS MANAGEABLE --> Is this even OK to do? probably yes just to get weightings
            # print(FE)
            ene_spline = make_interp_spline(betas, energy, k=3)
            ene_sq_spline = make_interp_spline(betas, energy_sq, k=3)
            FE_spline = make_interp_spline(betas, FE, k=3)
            ene_interp = ene_spline(beta_model)
            ene_sq_interp = ene_sq_spline(beta_model)
            ene_matrix[i,:] = ene_interp
            ene_sq_matrix[i,:] = ene_sq_interp
            obs_matrix[i,:] = (ene_sq_interp - (ene_interp*L)**2)*beta_model**2
            FE_interp[i,:] = FE_spline(beta_model)

            for j in range(N_steps):
                FE_weights[i,j] = np.exp(FE_interp[i,j] - max(FE_interp[:,j]))
                FE_weight_denoms[j] += FE_weights[i,j]
    
        for i in range(N_steps):
            FE_weights[:,i] /= FE_weight_denoms[i]
        
        unbiased_obs = np.zeros(N_steps)
        unbiased_ene = np.zeros(N_steps)
        unbiased_ene_sq = np.zeros(N_steps)
        for i in range(len(files)):
            unbiased_ene += ene_matrix[i,:] * FE_weights[i,:]
            unbiased_ene_sq += ene_sq_matrix[i,:] * FE_weights[i,:]
            
        unbiased_obs = (unbiased_ene_sq - (unbiased_ene*L)**2)*beta_model**2
        

    ### Bootstrap method
    B = 100

    
    bootstrap_obs = np.zeros(shape = (B,N_steps))
    se_obs_bootstrap = np.zeros(shape = N_steps)
    bootmean_std = np.zeros(len(unbiased_obs))
    for i in range(N_steps):
        sample = obs_matrix[:,i]
        boot_means = []
        for _ in range(B):
            bootsample = np.random.choice(sample,size=len(sample), replace=True)
            boot_means.append(bootsample.mean())
        bootmean_std[i] = np.std(boot_means,ddof = 1)
    se_obs_bootstrap = bootmean_std


    # for i in range(B):
    #     for j in range(N_steps):
    #         for k in range(N_samples):
    #             r = np.random.randint(0,len(files))
    #             bootstrap_obs[i,j] += obs_matrix[r,j]
    # bootstrap_obs /= N_samples

    # # print(bootstrap_obs)
    # for i in range(N_steps):
    #     se_obs_bootstrap[i] = np.sqrt((np.std(bootstrap_obs[:,i])**2)/(B-1))
    # print(B)

    # for i in range(B):
    #     for j in range(N_steps):
    #         for k in range(N_samples):
    #             r = np.random.randint(0,len(files))
    #             bootstrap_obs[i,j] += obs_matrix[r,j]
    # bootstrap_obs /= N_samples



    return unbiased_obs, se_obs_bootstrap



# getUnbiasedEstimate(kappa=0.6,L=32,R=1000,mode="t",quantity=6)

def getUnbiasedFE(kappa,L,R,mode):
    kappastr = f"{kappa:.2f}"
    modestr = "Two-Replica_Method" if (mode == "t") else "Wolff_Method"
    dirpath = f"/Users/shanekeiser/Downloads/production-run/{modestr}/{kappastr}_kappa/{L}_L/{R}_R"
    filepattern = f"emcx_data_*"
    files = glob.glob(f"{dirpath}/{filepattern}")
    beta_model = np.linspace(0.1,2.5,N_steps)   
    FE_interp_vals = np.empty(shape = (len(files), N_steps))
    unbiased_FE = np.empty(N_steps)
    for i in range(len(files)):
        df = pd.read_csv(files[i])
        betas = df["Beta"]
        FE = df["Free Energy"] # need to manipulate F a little, based on what is calculated in C++
        FE_interp = make_interp_spline(betas,FE,k=3)
        FE_model = FE_interp(beta_model)
        FE_interp_vals[i,:] = FE_model
    
    for i in range(N_steps):
        unbiased_FE[i] = np.log((1/len(files))*np.sum(np.exp(FE_interp_vals[:,i])))
        unbiased_FE[i] *= -1/beta_model[i]
    return unbiased_FE


def getProbFindGS(kappa,L,R,quantity):
    kappastr = f"{kappa:.2f}"
    dirpath_t = f"/Users/shanekeiser/Downloads/production-run/Two-Replica_Method/{kappastr}_kappa/{L}_L/{R}_R"
    dirpath_w = f"/Users/shanekeiser/Downloads/production-run/Wolff_Method/{kappastr}_kappa/{L}_L/{R}_R"
    filepattern = f"emcx_data_*"
    filepattern_2 = f"param_info_*"
    files_t = glob.glob(f"{dirpath_t}/{filepattern}")
    files_w = glob.glob(f"{dirpath_w}/{filepattern}")
    df_t = pd.read_csv(files_t[0])
    df_w = pd.read_csv(files_w[0])
    # p_info = np.loadtxt(glob.glob(f"{dirpath}/{filepattern_2}")[0])
    header = list(df_t)[quantity]

    final_val = 0 # Placeholder

    if quantity == 11 and kappa > 0.5: # Which other quantities can we use?
        final_val = 0.25

    gs_found_counter_t = 0
    gs_found_counter_w = 0
    obs_array_t = np.empty(shape = (len(files)))
    obs_array_w = np.empty(shape = (len(files)))

    for i in range(len(files)):
        obs_array_t[i] = getObs(files_t[i],quantity)[-1]
        if np.abs(obs_array_t[i] - final_val) < 1e-4: # 1e-4 tolerance
            gs_found_counter_t += 1
        obs_array_w[i] = getObs(files_w[i],quantity)[-1]
        if np.abs(obs_array_w[i] - final_val) < 1e-4: # 1e-4 tolerance
            gs_found_counter_w += 1
    
    prob_t = gs_found_counter_t/len(files)
    prob_w = gs_found_counter_w/len(files)
    
    print("Probability to reach f = 0.25 (Dominant frequency)")
    print(f"R = {R}, L = {L}, kappa = {kappa}")
    print(f"For TR: {prob_t} | For Wolff: {prob_w}")


# comparer = 11
# compareQuantity(kappa = 0.6, L = 32, R = 1000, quantity = comparer)
# compareQuantity(kappa = 0.6, L = 32, R = 2000, quantity = comparer)
# compareQuantity(kappa = 0.6, L = 32, R = 5000, quantity = comparer)

# getProbFindGS(kappa=0.6,L=32,R=1000,quantity=11)
# getProbFindGS(kappa=0.6,L=32,R=2000,quantity=11)
# getProbFindGS(kappa=0.6,L=32,R=5000,quantity=11)

kappa = 0.6
L = 48
R = 1600
quantity = 6
mode = "t"
modestr = "Two Replica Method" if (mode == "t") else "Wolff Method"
modestr = "Two-Replica_Method" if (mode == "t") else "Wolff_Method"
dirpath = f"/Users/shanekeiser/Downloads/production-run/{modestr}/{kappastr}_kappa/{L}_L/{R}_R"
filepattern = f"emcx_data_*"
files = glob.glob(f"{dirpath}/{filepattern}")
# print(files)
df_example = pd.read_csv(files[0])
header = list(df_example)[quantity]
print(header)


R_vals = np.array([80,800,1600,16000]) # R*S = 160000, epsilon = 0.1

unbiased_observables = np.empty(shape = (N_steps,len(R_vals)))
standard_errors = np.empty(shape = (N_steps,len(R_vals)))
for i in range(len(R_vals)):
    unbiased_observables[:,i], standard_errors[:,i] = getUnbiasedEstimate(kappa = kappa, L = L, R = R_vals[i], mode = mode, quantity = quantity)

# Customize plot appearance
plt.rc('text', usetex=False)  # Use LaTeX for all text
plt.rc('font', family='serif', size=10)  # Use serif fonts with size 10
plt.rc('axes', labelsize=12)  # Axis label size
plt.rc('legend', fontsize=10)  # Legend font size

fig, ax1 = plt.subplots(figsize=(9, 6))  # Adjust the figure size for the journal aspect ratio

from matplotlib.cm import get_cmap
# colors = get_cmap("Set1", len(R_vals)+1)
colors =  np.array(["#fa0707", "#363636", "#30fc03", "#0207fa","#f757f2","#0dffef"])#, "#f757f2", "#0dffef", "#ffaa0d"])
ecolors = np.array(["#960404", "#050505", "#1c9701", "#01037d", "#f757e3","#0dff13"])
markers = ["o", "v", "^", "s", "P", "X"]

inset = False
if inset == True:
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    inset_ax = inset_axes(ax1, width="40%", height="50%", loc='upper left')
    inset_ax.set_xlim(0.92, 1.21)  # Zoom in on the region of interest
    # inset_ax.set_xlim(0.76, 1.06)  # Zoom in on the region of interest
    inset_ax.set_ylim(9.9/48, 12.1/48)
    # inset_ax.set_ylim(8.9/48, 10.1/48)
    inset_ax.grid(visible=True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
    inset_ax.tick_params(left=False, right=False, top=False, bottom=False, labelleft=False, labelbottom=False)

# S_times_R = 160000
# M_vals = np.array([194,155,191,184,187,196])
M_vals = np.array([199,199,198,196])
# M_vals = np.array([190,194,191,184,190])
# S_vals = np.array([200,100,50,20,10])
S_vals = 160000/R_vals
e_vals = np.full(len(R_vals),0.05)
# e_vals = np.array([0.1,0.1,0.05,0.15,0.1])
# Plot the first dataset with error bars
for i in range(len(R_vals)):
    label = fr"R = {R_vals[i]}, S = {S_vals[i]}, M = {M_vals[i]}, $\epsilon$ = {e_vals[i]}"
    ax1.errorbar(beta_model, unbiased_observables[:,i], standard_errors[:,i], 
                linewidth=0.1, elinewidth=0.8, 
                marker=markers[i], markersize=1.5, 
                color='none', ecolor=ecolors[i], 
                mec = colors[i], mew = 0.5,
                label=label)
    if inset == True:
        inset_ax.errorbar(beta_model, unbiased_observables[:,i], standard_errors[:,i], 
                linewidth=0.1, elinewidth=0.8, 
                marker=markers[i], markersize=1, 
                color='none', ecolor=ecolors[i], 
                mec = colors[i], mew = 0.5,
                label=label)


# Axis labels
ax1.set_xlabel(r"Inverse temperature $\beta$", fontsize=12)
ax1.set_ylabel(f"{header}", fontsize=12)

# Limits and grid
ax1.set_xlim(0.2, 1.4)
ax1.set_xlim(0.9, 1.3)
ax1.set_ylim(0,3)
if quantity == 11:
    ax1.set_ylim(0.20,0.26)

ax1.grid(visible=True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)

# Title with LaTeX formatting
ax1.set_title(rf"Unbiased Estimator with Standard Error for {header}" +f"\n"
              + rf"$\kappa = {kappa}$, $L = {L}$", fontsize=11)

# Legend
ax1.legend(loc='upper left', frameon=True, markerscale = 4)

from fractions import Fraction
if quantity == 11:
    ax1.hlines(y = np.arange(((L/8)+1)/L, ((L/4)+1)/L, 1/L), xmin = 0, xmax = 2.5, linewidth = 1, linestyle = 'dashed', alpha = 0.6, color = "black")
    # ax1.set_title(rf"Unbiased Estimator with Standard Error for Dominant Wavenumber" +f"\n"
    #           + rf"$\kappa = {kappa}$, $L = {L}$, {modestr}", fontsize=11)
    ax1.set_ylabel(r"Wavenumber$/2\pi$", fontsize=12)
    if inset == True:
        inset_ax.hlines(y = np.arange(((L/8)+1)/L, ((L/4)+1)/L, 1/L), xmin = 0, xmax = 2.5, linewidth = 1, linestyle = 'dashed', alpha = 0.6, color = "black")
        ax1.set_yticks(ticks = np.arange(((L/8)+1)/L, ((L/4)+1)/L, 1/L), labels = (f'{str(int(i*L))}/{str(L)}' for i in np.arange(((L/8)+1)/L, ((L/4)+1)/L, 1/L)))
# Optimize layout and save


# # plt.tight_layout()
# from datetime import date
# datestr = date.today().strftime("%d%b%y")
# descriptor = f"{header}-{modestr}"
# print(datestr)
# if kappa == 0.6:
#     kappafilestr = "k-0-6"
# elif kappa == 0:
#     kappafilest = "k-0"
# plt.savefig(f"/Users/shanekeiser/Documents/ANNNI/figures/{datestr}_{descriptor}_{kappafilestr}_L-{L}_R-comparison.png", dpi=300)  # Save in high resolution
plt.show()
