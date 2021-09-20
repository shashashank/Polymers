import numpy as np
import random
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from numba import prange,jit
import time

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

# Evolves the points of the polymer due to thermal fluctuation
@jit(fastmath=True, nopython=True)
def animate(pos, phi, chain_length, bond_length=2, time_step=1e-3, std_deviation=1):
    n = chain_length        # Length of the polymer
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



@jit(fastmath=True, parallel = True)
def yes(array1, array2, time, step, length):
    L = array1
    L1 = array2
    x = np.zeros(time)
    for i in prange(L.shape[0]):
        n = 10 + i
        pos = np.zeros((n,2))   # position array
        phi = np.zeros(1)     # angle

        # Setting the IC
        for f in range(len(pos)):
            pos[f,0] = f * length
        
        for j in range(time):
            x[j] = animate(pos, phi, n, length, step, 1)
            if j==0: L1[i] = x[j]
        L[i] = np.sum(x)/time
    return L
        #f = int(np.floor(n/2))
        #box = 20
        # plt.axis((pos[f,0] - box, pos[f,0] + box, pos[f,1] - box, pos[f,1]+ box))
        # for k in range(n):
        #     plt.plot(pos[k,0], pos[k,1], 'bo', markersize=5)
        #     if k > 0: plt.plot([pos[k-1,0],pos[k,0]], [pos[k-1,1],pos[k,1]], 'gray', linestyle=':', marker='')
        # plt.show()

n = 100
l = 2                  # bond length
t = int(1e8)           # simulation time
s = 1e-5               # step length
L = np.zeros(n)        
L1 = np.zeros(n)       # shows initial E2E dist

start = time.time()
y = yes(L, L1, t, s, l)
x = np.linspace(10, len(y)+10, len(y))
plt.plot(x, y, 'b')
plt.plot(x, L1, 'r')
plt.tight_layout()
plt.show()
end = time.time()
print("Elapsed (after compilation) = %s" % (end - start))