# ADD LATER 

# Import libraries 
import numpy as np 
import matplotlib.pyplot as plt 
import time 

import matplotlib.animation as anim 
from matplotlib.animation import FuncAnimation

start_time = time.time()

# Parameters 
N = 100 

# Load files 
eta_1 = np.load("results/xi_1.npy")[::3, :]
time_list = np.linspace(0,200000)

fig, axs = plt.subplots(1,1)

def animate(i):

    axs.clear()

    axs.plot(range(N), eta_1[i,:])
    axs.set_ylim([0,0.001])

    print(i)


ani = FuncAnimation(fig, animate, frames=3000, interval=1, repeat=False, save_count=3000)
Writer = anim.writers['ffmpeg']
writer = Writer(fps=60, bitrate=2000)

ani.save(filename="results/animation.mp4", writer=writer, dpi=300)
