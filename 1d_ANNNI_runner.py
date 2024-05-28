"""   
    This file was used to remotely run an early 1D simulation, before I figured out
    you could very easily (and more efficiently) implement Temperature and kappa
    iteration directly in the C++ file.

    But this was a good practice in using subprocess and os.

    I leave notes at the bottom which indicate the next steps for this idea.
"""


import numpy as np
import subprocess
import os
import signal
import matplotlib.pyplot as plt
from subprocess import Popen, PIPE
from time import sleep

T_vals = np.linspace(0,2,25)
k_vals = np.linspace(0,2,25)

program_path = "/Users/shanekeiser/Documents/Spring 2024/Machta Spring/cluster1d_take2_exe"

i = 0
for k in k_vals:
    i+= 1
    for T in T_vals:
        p = Popen([program_path], stdout=PIPE, stdin=PIPE)
        value1 = bytes(str(k) + '\n', 'UTF-8') 
        value2 = bytes(str(T) + '\n', 'UTF-8') 
        p.stdin.write(value1)
        p.stdin.write(value2)
        p.stdin.flush()
        sleep(0.1)
    print(f"{i} sets done!")

print('Completely done!')

### SO FAR, we have just done the part which runs the code, for any number of T and kappa values. 
### Next, we read in and analyze the data from each run of the code, and do a colormap plot for some quantity
### with T and kappa as the x and y axis.
        
### Note that one thing left to do is figure out how many steps is in one sweep, and whether or not the system
### is equilibriating yet.

### Rewrite the C++ file to only output the E,M,C,X average and std? and whatever other order parameter there is.

