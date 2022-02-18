# Compares eta and xi values to spin values to make sure that the normalisation is working 

import numpy as np 
import matplotlib.pyplot as plt 

spin_up_data = np.loadtxt("results/spin_up.dat")
spin_down_data = np.loadtxt("results/spin_down.dat")

eta_0_data = np.loadtxt("results/eta_0.dat")
xi_0_data = np.loadtxt("results/xi_0.dat")
eta_1_data = np.loadtxt("results/eta_1.dat")
xi_1_data = np.loadtxt("results/xi_1.dat")

time_list = spin_up_data[:,0]
spin_up_list = spin_up_data[:,1]
spin_down_list = spin_down_data[:,1]

eta_list = eta_0_data[:,1]
xi_list = xi_0_data[:,1]

for t in range(len(time_list)):

    eta_list[t] += sum(eta_1_data[t,:])
    xi_list[t] += sum(xi_1_data[t,:])

##### ----- Plotting -----

plt.plot(time_list, eta_list, label="eta", lw=2)
plt.plot(time_list, spin_down_list, label="spin down", ls="dotted", lw=2)

plt.legend()
plt.show()


