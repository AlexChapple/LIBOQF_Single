# Plots second order coefficients 

import numpy as np 
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt 
import matplotlib.cm as cm
import matplotlib

time_stamp_init_time = 5000 
time_stamp_duration = 5000 + 16000

eta_1 = np.loadtxt("results/eta_1.dat")
xi_1 = np.loadtxt("results/xi_1.dat")

time_list = np.linspace(0,5,200000)

shape = np.shape(eta_1)
x1 = []
y1 = []

dz = []
dz2 = []
for j in range(shape[1]):
    for k in range(shape[0]):
        if k % 100 == 0 and k > time_stamp_init_time and k < time_stamp_init_time + time_stamp_duration:
            x1.append(j)
            y1.append(time_list[k])
            dz.append(eta_1[k,j])
            dz2.append(xi_1[k,j])

print("sorted through data")

dx = 1 
dy = abs(time_list[21] - time_list[20])

z1 = np.zeros_like(x1)

fig1 = plt.figure(1)
axs1 = fig1.add_subplot(111, projection='3d')

cmap = cm.get_cmap('jet')
max_val = max(dz)
rgba = [cmap(i / max_val) for i in dz]

axs1.bar3d(x1,y1,z1,dx,dy,dz,shade=True,color=rgba)
axs1.set_xlabel("N")
axs1.set_ylabel("time")
axs1.set_zlabel("\n\n$\eta^{(1)}$")

fig2 = plt.figure(2)
axs2 = fig2.add_subplot(111, projection='3d')

cmap = cm.get_cmap('jet')
max_val = max(dz2)
rgba = [cmap(i / max_val) for i in dz2]

# axs2.bar3d(x1,y1,z1,dx,dy,dz2,shade=True,color=rgba)
axs2.plot_surface(x1, y1, dz, cmap=rgba, lw=0)
axs2.set_xlabel("N")
axs2.set_ylabel("time")
axs2.set_zlabel("\n\n$\\xi^{(1)}$")


plt.show()