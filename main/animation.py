import numpy as np
import random
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from numba import prange,jit

box = 20    #area of graph
n=15
pos = np.zeros([n,2])   #position array
phi = np.zeros([1])     #angle array
length = 2              #distance between the two particles
mks = 5                 #size of the two particles
pos[1,0] = 0            #arbitrary initial position
dt = 0.5                #delta_t
sqrt_dt = np.sqrt(dt)


# Generates the force based on the displacement between the two polymer atoms
@jit(fastmath=True)
def spring(x, y):
    k = 2
    distance = np.sqrt(x**2 + y**2)
    return -k*(distance - length)


def animate(i):
    plt.clf()
    plt.axis(( -1*box, box, -1*box, box))
    for i in range(n):
        # Stochastic Term
        pos[i,0] += np.random.normal(0, 0.5)*sqrt_dt
        pos[i,1] += np.random.normal(0, 0.5)*sqrt_dt

        # Spring Force Term
        if i!=(n-1):
            x = pos[i+1,0] - pos[i,0]
            y = pos[i+1,1] - pos[i,1]
            F = spring(x, y)
            phi = np.arctan2(y, x)
            pos[i,0] += F*np.cos(phi)*dt
            pos[i,1] += F*np.sin(phi)*dt
        if i!=0:
            x1 = pos[i,0] - pos[i-1,0]
            y1 = pos[i,1] - pos[i-1,1]
            F1 = spring(x1, y1)
            phi1 = np.arctan2(y1, x1)
            pos[i,0] += F1*np.cos(phi1)*dt
            pos[i,1] += F1*np.sin(phi1)*dt
        
        plt.plot(pos[i,0], pos[i,1], 'bo', markersize=mks)
        if i > 0:
            plt.plot([pos[i-1,0],pos[i,0]], [pos[i-1,1],pos[i,1]], 'gray', linestyle=':', marker='')

ani = FuncAnimation(plt.gcf(), animate, interval=100)
plt.tight_layout()
plt.show()
