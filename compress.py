# Takes in the .dat files from Fortran and turn them into npy files for quicker load times 

import numpy as np 

spin_up_data = np.loadtxt("results/spin_up.dat")
spin_down_data = np.loadtxt("results/spin_down.dat")
eta_0 = np.loadtxt("results/eta_0.dat")
xi_0 = np.loadtxt("results/xi_0.dat")
eta_1 = np.loadtxt("results/eta_1.dat")
xi_1 = np.loadtxt("results/xi_1.dat")

# Turn them into npy files 

np.save("results/spin_up.npy", spin_up_data)
np.save("results/spin_up.npy", spin_up_data)
np.save("results/eta_0.npy", eta_0)
np.save("results/xi_0.npy", xi_0)
np.save("results/eta_1.npy", eta_1)
np.save("results/xi_1.npy", xi_1)