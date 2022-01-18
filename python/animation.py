import numpy as np
import random
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from numba import prange,jit

box = 10                # side length of the plot
n=10                    # Length of the polymer
pos = np.zeros([n,2])   # position array
phi = np.zeros([1])     # angle
length = 1              # distance between the two particles
mks = 5                 # size of the two particles
pos[1,0] = 0            # arbitrary initial position
dt = 1e-2                # delta_t
sqrt_dt = np.sqrt(dt)


# Generates the f
# orce based on the displacement between the two polymer atoms
@jit(fastmath=True)
def spring(x, y):
    k = 20
    distance = np.sqrt(x**2 + y**2)
    return -k*(distance - length)

# Calculates the end-to-end distance of the polymer
@jit(fastmath=True)
def e2e_distance(p_init, p_fin):
    return np.sqrt((p_fin[0] - p_init[0])**2 + (p_fin[1] - p_init[1])**2)


def animate(i):
    plt.clf()
    plt.axis(( -box, box, -box, box))

    # Setting the IC
    j = - np.floor(len(pos)/2)
    for i in range(len(pos)):
        pos[i,0] = j*length
        j +=1

    for i in range(n):

        # Stochastic Term
        pos[i,0] += np.random.normal(0, 0.1)*sqrt_dt
        pos[i,1] += np.random.normal(0, 0.1)*sqrt_dt

        # Spring Force term
        if i!=(n-1):
            x = pos[i+1,0] - pos[i,0]
            y = pos[i+1,1] - pos[i,1]
        if i!=0:
            x = pos[i,0] - pos[i-1,0]
            y = pos[i,1] - pos[i-1,1]
        F = spring(x, y)
        phi = np.arctan2(y, x)      # Calculates the angle between two units
        pos[i,0] += F*np.cos(phi)*dt
        pos[i,1] += F*np.sin(phi)*dt

        #print("End to end distance is:%8.4f" % e2e_distance(pos[0,:], pos[n-1,:]))
        plt.plot(pos[i,0], pos[i,1], 'bo', markersize=mks)
        if i > 0:
            plt.plot([pos[i-1,0],pos[i,0]], [pos[i-1,1],pos[i,1]], 'gray', linestyle=':', marker='')

ani = FuncAnimation(plt.gcf(), animate, interval=1)
plt.tight_layout()
plt.show()
