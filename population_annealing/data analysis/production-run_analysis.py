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
def getUnbiasedEstimate(kappa,L,R,mode,quantity):
    kappastr = f"{kappa:.2f}"
    modestr = "Two-Replica_Method" if (mode == "t") else "Wolff_Method"
    dirpath = f"/Users/shanekeiser/Downloads/production-run/{modestr}/{kappastr}_kappa/{L}_L/{R}_R"
    filepattern = f"emcx_data_*"
    files = glob.glob(f"{dirpath}/{filepattern}")
    beta_model = np.linspace(0.1,2.5,N_steps)   
    FE_interp_vals = np.zeros(shape = (len(files), N_steps))
    FE_weights = np.zeros(shape = (len(files), N_steps))
    FE_weight_denom = np.zeros(shape = N_steps)
    #header = list(df)[quantity]
    for i in range(len(files)):
        df = pd.read_csv(files[i])
        betas = df["Beta"]
        FE = df["Free Energy"]
        FE_interp = make_interp_spline(betas,FE,k=3)
        FE_model = FE_interp(beta_model)
        FE_interp_vals[i,:] = FE_model
        FE_weights[i,:] = np.exp(FE_model)
        FE_weight_denom += np.exp(FE_model)
    
    for i in range(len(files)):
        FE_weights[i,:] /= FE_weight_denom
    print(FE_weights)
    observables = np.zeros(shape = (len(files), N_steps))
    unbiased_obs = np.zeros(shape = N_steps)
    for i in range(len(files)):
        observables[i,:] = getObs(files[i],quantity)
    
        unbiased_obs += FE_weights[i,:] * observables[i,:]

    # unbiased_obs = observables * FE_weights
    print(unbiased_obs)
    plt.plot(beta_model, unbiased_obs)
    plt.plot(beta_model, FE_model)
    plt.show()

getUnbiasedEstimate(kappa=0.6,L=32,R=1000,mode="t",quantity=6)

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


comparer = 11
compareQuantity(kappa = 0.6, L = 32, R = 1000, quantity = comparer)
compareQuantity(kappa = 0.6, L = 32, R = 2000, quantity = comparer)
compareQuantity(kappa = 0.6, L = 32, R = 5000, quantity = comparer)

getProbFindGS(kappa=0.6,L=32,R=1000,quantity=11)
getProbFindGS(kappa=0.6,L=32,R=2000,quantity=11)
getProbFindGS(kappa=0.6,L=32,R=5000,quantity=11)


# Unbiased Estimator of FE

unbiased_FE_t = getUnbiasedFE(kappa = 0.6, L = 32, R = 1000, mode = "t")
unbiased_FE_p = getUnbiasedFE(kappa = 0.6, L = 32, R = 1000, mode = "p")

plt.plot(beta_model, unbiased_FE_t, label = "Two Replica")
plt.plot(beta_model, unbiased_FE_p, label = "Wolff")
plt.title("Comparison: Unbiased estimator of Free Energy")
plt.legend()
plt.tight_layout()
plt.show()

# Compare FE over the two methods and R
R_vals = np.array([1000,2000,5000])

t_FE_stats = np.empty(shape = (len(R_vals), 2))
p_FE_stats = np.empty(shape = (len(R_vals), 2))

for i in range(len(R_vals)):
    t_FE_stats[i,:] = getMeanFE(kappa = 0.6, L = 32, R = R_vals[i], mode = "t")
    p_FE_stats[i,:] = getMeanFE(kappa = 0.6, L = 32, R = R_vals[i], mode = "p")

plt.errorbar(R_vals, t_FE_stats[:,0], t_FE_stats[:,1], marker = '.', label = "TR", elinewidth = 0.5, linewidth = 0, capsize = 0.2)
plt.errorbar(R_vals, p_FE_stats[:,0], p_FE_stats[:,1], marker = '.', label = "Wolff", elinewidth = 0.5, linewidth = 0, capsize = 0.2)
plt.legend()
plt.title(r"Comparison of Free energy per spin $\tilde{F}_{1}/N = - \frac{1}{N\beta_1} \sum_{l=K}^{1} \ln Q(\beta_l, \beta_{l-1})$," + f"\n" +  r"interpolated at $\beta \approxeq 2$")
plt.xlabel(r"Population size $R$")
plt.ylabel(r"$\tilde{F}_{1}/N$")
plt.tight_layout()
plt.show()
