import numpy as np
import matplotlib.pyplot as plt
import pandas as pd



def importData(kappastr, size, index):
    fname = f"/Users/shanekeiser/Documents/ANNNI/populationannealing/data/t/emcx_data_{kappastr:.2f}_kappa_{size}_L.csv"
    info_name = f"/Users/shanekeiser/Documents/ANNNI/populationannealing/data/t/parameter_info_{kappastr:.2f}_kappa_{size}_L.csv"
    ene_h_name = f"/Users/shanekeiser/Documents/ANNNI/populationannealing/data/t/ene_hist_{kappastr:.2f}_kappa_{size}_L.csv"
    mag_h_name = f"/Users/shanekeiser/Documents/ANNNI/populationannealing/data/t/mag_hist_{kappastr:.2f}_kappa_{size}_L.csv"

    cols = np.arange(1,3000,1)
    info = np.loadtxt(info_name, delimiter = ',', dtype = str)
    df = pd.read_csv(fname)
    ene_hist = np.loadtxt(ene_h_name, delimiter = ',', dtype = float, usecols = cols)
    mag_hist = np.loadtxt(mag_h_name, delimiter = ',', dtype = float, usecols = cols)
    betas = df["Beta"]
    ene_slice = ene_hist[index]
    ene_plotter = ene_slice[ene_slice < 5000]
    mag_slice = mag_hist[index]
    mag_plotter = mag_slice[mag_slice < 5000]
    fig, axs = plt.subplots(2,1)
    
    print(sorted(ene_plotter)[1])
    print(sorted(mag_plotter)[1])
    print(betas[index])

    axs[0].hist(ene_plotter, bins = 30)
    axs[0].set_title(f"Energy histogram for kappa = {kappastr}, L = {size}, beta = {betas[index]}")

    axs[1].hist(mag_plotter, bins = 30)
    axs[1].set_title(f"Magnetization histogram for kappa = {kappastr}, L = {size}, beta = {betas[index]}")

    plt.tight_layout()
    plt.show()


importData(0,32,117) # Here, index = 131 corresponds to the most devious spike

