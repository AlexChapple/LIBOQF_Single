import numpy as np 
import matplotlib.pyplot as plt 

spin_up_data = np.loadtxt("results_julia/spin_up.dat")

time_list = spin_up_data[:,0]
spin_up_list = spin_up_data[:,1]

plt.plot(time_list, spin_up_list, lw=2)

plt.show()