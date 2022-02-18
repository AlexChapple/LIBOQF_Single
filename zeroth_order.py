# Plot first order coefficients 

import numpy as np 
import matplotlib.pyplot as plt 

eta_0_data = np.load("results/eta_0.npy")
xi_0_data = np.load("results/xi_0.npy")

time_list = eta_0_data[:,0]
eta_0_list = eta_0_data[:,1]
xi_0_list = xi_0_data[:,1]

plt.plot(time_list, eta_0_list, label="$\eta^{(0)}$", lw=2)
plt.plot(time_list, xi_0_list, label="$\\xi^{(0)}$", lw=2)

plt.legend()
plt.show()