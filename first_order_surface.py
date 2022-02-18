# Plots first order coefficients 


import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import cm 
from mpl_toolkits.mplot3d import Axes3D
import matplotlib

N = 100 
time_stamp_init_time = 37000 
time_stamp_end_time = time_stamp_init_time + 20000

eta_1 = np.load("results/eta_1.npy")[time_stamp_init_time:time_stamp_end_time, :]
xi_1 = np.load("results/xi_1.npy")[time_stamp_init_time:time_stamp_end_time, :]

time_list = np.linspace(0,5,200000)[time_stamp_init_time:time_stamp_end_time]
N_range = range(N)

shape = np.shape(eta_1)
x1 = []
y1 = []

z1 = []
z2 = []
for j in range(shape[1]):
    for k in range(shape[0]):
        x1.append(j)
        y1.append(time_list[k])
        z1.append(eta_1[k,j])
        z2.append(xi_1[k,j])

print("sorted through data")

### Reorganise arrays 

X = np.reshape(x1, (N, time_stamp_end_time - time_stamp_init_time))
Y = np.reshape(y1, (N, time_stamp_end_time - time_stamp_init_time))
Z1 = np.reshape(z1, (N, time_stamp_end_time - time_stamp_init_time))
Z2 = np.reshape(z2, (N, time_stamp_end_time - time_stamp_init_time))


### Plotting 
fig = plt.figure(1)
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, Z1, cmap=cm.jet)
ax.set_xlabel("N")
ax.set_ylabel("$\gamma t$")
ax.set_zlabel("\n\n$\eta^{(1)}$")

fig = plt.figure(2)
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, Z2, cmap=cm.jet)
ax.set_xlabel("N")
ax.set_ylabel("$\gamma t$")
ax.set_zlabel("\n\n$\\xi^{(1)}$")

Z3 = Z1 + Z2
fig = plt.figure(3)
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, Z3, cmap=cm.jet)
ax.set_xlabel("N")
ax.set_ylabel("$\gamma t$")
ax.set_zlabel("\n\n$\eta^{(1)} + \\xi^{(1)}$")

plt.show()