import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import glob

relative_path = "../data/t/"
### Single Kappa/size measure

kappa = 0
def Analyze(kappa = None, size = None, mode = 't', quantities = [], normalize = False, marker = '.'):
    while kappa == None or size == None:
        print("Kappa or size missing. Please input both. ")
        kappa = float(input("Kappa: "))
        size = input("Length (16,32,64,128,256): ")
    kappastr = f'{kappa:.2f}'
    df = f"/Users/shanekeiser/Documents/ANNNI/populationannealing/data/" + mode + f"/emcx_data_{kappastr}_kappa_{size}_L.csv"
    param_info = f"/Users/shanekeiser/Documents/ANNNI/populationannealing/data/" + mode + f"/parameter_info_{kappastr}_kappa_{size}_L.csv"

    makePlots(df, param_info, mode, quantities, normalize = normalize, marker = marker)
    return 0

# Quantities (in order):
# 0: Beta    1: Energy   2: Energy Squared   3: Mag   4: Mag Squared   5: Absolute Mag  6: Specific Heat
# 7: Susceptibility    8: Cluster Size   9: Non-Wrapping Cluster Size   10: Wrapping Probability

quantities = { 0 : "Beta",
               1 : "Energy",
               2 : "Energy Squared",
               3 : "Magnetization",
               4 : "Magnetization Squared",
               5 : "Absolute Magnetization",
               6 : "Specific Heat",
               7 : "Susceptibility",
               8 : "Cluster Size",
               9 : "Non-Wrapping Cluster Size",
               10: "Wrapping Probability",
               11: "Dominant Frequency",
               12: "Dominant Amplitude",
               13: "Rho T",
               14: "No. of Families",
               15: "Overlap",
               16: "Absolute Overlap",
               17: "Overlap Variance",
               }


def makePlots(fname = "nil", info_name = "nil", mode = "t", quantities = [], normalize = False, marker = '.'):
    if fname == "nil" or info_name == "nil":
        print("Please include both the data and parameter files. Exiting...")
        return (1)
    info = np.loadtxt(info_name, delimiter = ',', dtype = str, skiprows=1)
    df = pd.read_csv(fname)
    
    NN = float(info[1])**2
    if normalize == True:
        df["Cluster Size"] = df["Cluster Size"].div(NN)
        df["Energy Squared"] = df["Energy Squared"].div(NN)
        df["Magnetization Squared"] = df["Magnetization Squared"].div(NN)
        df["Non-Wrapping Cluster Size"] = df["Non-Wrapping Cluster Size"].div(NN)
        df["Susceptibility"] = df["Susceptibility"].div(NN)
        df["Overlap Variance"] = df["Overlap Variance"].div(NN) 
        
    if sum(np.square(quantities)) == 130:
        NWCS = np.array(df["Non-Wrapping Cluster Size"])
        B = np.array(df["Beta"])
        new_NWCS = NWCS*B
        df["Non-Wrapping Cluster Size"]=new_NWCS
    
    modestr = "nil"
    if mode == "t":
        modestr = "Two replica"
    elif mode == "p":
        modestr = "Wolff"
    titlestr = f"kappa = {info[0]}, L = {info[1]}, INIT_POP_SIZE = {info[2]}, CULLING_FRAC = {info[3]}, algorithm = {modestr}"
    if quantities == []:

        title = f"Quantities of interest for\n{titlestr}"
        axes = df.plot(x='Beta', subplots = True, layout = [6,4], figsize = [10,8], title = title, legend = False, marker = marker, linewidth = 0)
        ax_x, ax_y = 0,0
        for i in range(1,len(list(df))):
            axes[ax_x, ax_y].set_ylabel(list(df)[i], fontsize = 'small')
            axes[ax_x, ax_y].grid()
            ax_y += 1
            if ax_y == 4:
                ax_x += 1
                ax_y = 0
        plt.tight_layout()
        plt.show()
    else:
        headers = []
        for i in quantities:
            headers.append(list(df)[i])
        title = f"Comparison plot for\n{titlestr}"
        df.plot(x = 'Beta', y = headers, title = title, marker = marker, linewidth = 0)
        plt.show()
    print(f"No. of Temperature steps = {len(df['Beta'])}")
    print(df.head(5))




# To compare across several kappa or lengths
def Compare(kappas = [], sizes = [], quantity = 0, normalize = False, marker = '.'):
    dfs = []
    infos = []

    while quantity == 0:
        print("No quantity provided. Please use keyword argument 'quantity', as such:\n \
               1 : Energy, \n2 : Energy Squared, \n3 : Magnetization, \n \
               4 : Magnetization Squared, \n5 : Absolute Magnetization, \n \
               6 : Specific Heat, \n7 : Susceptibility, \n \
               8 : Cluster Size, \n9 : Non-Wrapping Cluster Size, \n \
               10: Wrapping Probability")
        quantity = int(input("Desired quantity: "))

    kappastrs = []
    for kappa in kappas:
        kappastrs.append(f'{kappa:.2f}')

    #### COLOR MAPPING
    colors = plt.cm.turbo(np.linspace(0,1,np.max(a = (len(kappas), len(sizes)))))

    if len(sizes) == 0 and len(kappas) == 0:
        print("Please provide kappa or size value to begin analysis.")
        return 2
    elif len(sizes) == 1:
        for i in range(len(kappas)):
            fname = f"/Users/shanekeiser/Downloads/data/emcx_data_{kappastrs[0]}_kappa_{sizes[i]}_L.csv"
            info_name = f"/Users/shanekeiser/Downloads/data/parameter_info_{kappastrs[0]}_kappa_{sizes[i]}_L.csv"
            infos.append(np.loadtxt(info_name, delimiter = ',', dtype = str))
            dfs.append(pd.read_csv(fname))
        header = list(dfs[0])[quantity]
        print(header)
        for i in range(len(kappas)):
            df = dfs[i]
            info = infos[i]
            if normalize == True:
                NN = float(info[1])**2
                df["Cluster Size"] = df["Cluster Size"].div(NN)
                df["Energy Squared"] = df["Energy Squared"].div(NN)
                df["Magnetization Squared"] = df["Magnetization Squared"].div(NN)
                df["Non-Wrapping Cluster Size"] = df["Non-Wrapping Cluster Size"].div(NN)
                df["Susceptibility"] = df["Susceptibility"].div(NN)

            
            # Make plots
            if i == 0:
                titlestr = f'{header} comparison over varying $\kappa$ for $L$ = {info[1]}\nINIT_POP_SIZE = {info[2]}, CULLING_FRAC = {info[3]}'
                if normalize == True:
                    titlestr = titlestr + '\n(values normalized: divided by no. of spins)'
                ax = df.plot(x="Beta", y = header, label = f'$\kappa$ = {info[0]}', title = titlestr, color = colors[0], marker = marker)
            else:
                df.plot(ax = ax, x = "Beta", y = header, label = f'$\kappa$ = {info[0]}', color = colors[i], marker = marker)
        
        plt.show()
    elif len(kappas) == 1:
        for i in range(len(sizes)):
            fname = f"/Users/shanekeiser/Downloads/data/emcx_data_{kappastrs[0]}_kappa_{sizes[i]}_L.csv"
            info_name = f"/Users/shanekeiser/Downloads/data/parameter_info_{kappastrs[0]}_kappa_{sizes[i]}_L.csv"
            infos.append(np.loadtxt(info_name, delimiter = ',', dtype = str))
            dfs.append(pd.read_csv(fname))
        header = list(dfs[0])[quantity]
        # print(header)
        for i in range(len(sizes)):
            df = dfs[i]
            info = infos[i]
            if normalize == True:
                NN = float(info[1])**2
                df["Cluster Size"] = df["Cluster Size"].div(NN)
                df["Energy Squared"] = df["Energy Squared"].div(NN)
                df["Magnetization Squared"] = df["Magnetization Squared"].div(NN)
                df["Non-Wrapping Cluster Size"] = df["Non-Wrapping Cluster Size"].div(NN)
                df["Susceptibility"] = df["Susceptibility"].div(NN)
            # Make plots    
            if i == 0:
                titlestr = f'{header} comparison over varying $L$ for $\kappa$ = {info[0]}\nINIT_POP_SIZE = {info[2]}, CULLING_FRAC = {info[3]}'
                if normalize == True:
                    titlestr = titlestr + '\n(values normalized: divided by no. of spins)'
                ax = df.plot(x="Beta", y = header, label = f'$LEN$ = {info[1]}', title = titlestr, color=colors[0], marker = marker)
            else:
                df.plot(ax = ax, x = "Beta", y = header, label = f'$LEN$ = {info[1]}', color = colors[i], marker = marker)
        plt.show()

# Analyze(0.25,16, quantities = [8,9])
# Analyze(0.25,32, quantities = [8,9])

Analyze(0,16, mode = 'p', quantities = [], normalize=False)
# Analyze(0.6, 16, quantities = [15,16,17], normalize=True)    

# Analyze(0.25,64, normalize = True, quantities = [7,10])


# Compare(kappas = [0,0.25,0.5,0.75,1], sizes = [32], quantity = 7)


# Compare(kappas = [0.25], sizes = [16,32], quantity = 1, normalize = True)
