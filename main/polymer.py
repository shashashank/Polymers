import numpy as np
import random
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from numba import prange,jit

# Generates the force based on the displacement between the two polymer atoms
@jit(fastmath=True, nopython=True)
def spring(x, y, length):
    k = 20
    distance = np.sqrt(x*x + y*y)
    return -k*(distance - length)

# Calculates the end-to-end distance of the polymer
@jit(fastmath=True, nopython=True)
def e2e_distance(p_init, p_fin):
    return np.sqrt((p_fin[0] - p_init[0])*(p_fin[0] - p_init[0]) + (p_fin[1] - p_init[1])*(p_fin[1] - p_init[1]))

@jit(fastmath=True, nopython=True)
def animate(pos, phi, chain_length, bond_length=2, time_step=1e-3, std_deviation=1):
    n=chain_length          # Length of the polymer
    length = bond_length    # distance between the two particles
    dt = time_step          # delta_t
    sigma = std_deviation   # std dev of noise term
    sqrt_dt = np.sqrt(dt)

    for i in range(n):
        # Stochastic Term
        pos[i,0] = pos[i,0] + np.random.normal(0, sigma)*sqrt_dt
        pos[i,1] = pos[i,1] + np.random.normal(0, sigma)*sqrt_dt

        # Spring Force term
        if i!=(n-1):
            x = pos[i+1,0] - pos[i,0]
            y = pos[i+1,1] - pos[i,1]
        if i!=0:
            x = pos[i,0] - pos[i-1,0]
            y = pos[i,1] - pos[i-1,1]
        F = spring(x, y, length)
        phi = np.arctan2(y, x)      # Calculates the angle between two units
        pos[i,0] = pos[i,0] + F*np.cos(phi)*dt
        pos[i,1] = pos[i,1] + F*np.sin(phi)*dt

    return e2e_distance(pos[0,:], pos[n-1,:])

L = np.zeros(1000)
#@jit(fastmath=True, nopython=True, parallel = True)
def yes():
    length = 2                  # bond length
    time = int(1e5)                  # simulation time
    x = np.zeros(time)
    for i in range(len(L)):
        n = 10 + i
        pos = np.zeros([n,2])   # position array
        phi = np.zeros([1])     # angle

        # Setting the IC
        for f in range(len(pos)):
            pos[f,0] = f * length
        for j in range(time):
            x[j] = animate(pos, phi, n, length, 1e-3, 1)
        L[i] = np.sum(x)/time
        #f = int(np.floor(n/2))
        #box = 20
        # plt.axis((pos[f,0] - box, pos[f,0] + box, pos[f,1] - box, pos[f,1]+ box))
        # for k in range(n):
        #     plt.plot(pos[k,0], pos[k,1], 'bo', markersize=5)
        #     if k > 0: plt.plot([pos[k-1,0],pos[k,0]], [pos[k-1,1],pos[k,1]], 'gray', linestyle=':', marker='')
        # plt.show()

yes()
plt.plot(L)
plt.show

plt.tight_layout()
plt.show()
