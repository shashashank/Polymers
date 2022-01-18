import time
import random
import numpy as np
import matplotlib.pyplot as plt
from numba import prange,jit

# Generates the force based on the displacement between the two polymer atoms
@jit(fastmath=True, nopython=True)
def spring(X, Y, length, K=500):
    distance = np.sqrt(X*X + Y*Y)
    return -K*(distance - length)

# Calculates the end-to-end distance of the polymer
@jit(fastmath=True, nopython=True)
def e2e_distance(p_init, p_fin):
    X = (p_fin[0] - p_init[0])*(p_fin[0] - p_init[0])
    Y = (p_fin[1] - p_init[1])*(p_fin[1] - p_init[1])
    return np.sqrt(X + Y)

# Evolves the points of the polymer due to thermal fluctuation
@jit(fastmath=True, nopython=True, parallel = True)
def animate(pos, chain_length, bond_length=1, time_step=1e-3, std_deviation=1, K=500):
    n = chain_length        # Length of the polymer
    length = bond_length    # distance between the two particles
    sigma = std_deviation   # std dev of noise term
    k = K                   # spring constant
    dt = time_step          # delta_t
    sqrt_dt = np.sqrt(dt)
    
    # Stochastic Term
    for i in prange(n):
        pos[i,0] = pos[i,0] + np.random.normal(0.0, sigma)*sqrt_dt
        pos[i,1] = pos[i,1] + np.random.normal(0.0, sigma)*sqrt_dt
        
    # Spring Force term
    for i in range(n):
        F1 = F2 = 0
        if i!=(n-1):
            x1 = pos[i+1,0] - pos[i,0]
            y1 = pos[i+1,1] - pos[i,1]
            F1 = spring(x1, y1, length, k)
            phi1 = np.arctan2(y1, x1)
        if i!=0:
            x2 = pos[i-1,0] - pos[i,0]
            y2 = pos[i-1,1] - pos[i,1]
            F2 = spring(x2, y2, length, k)
            phi2 = np.arctan2(y2, x2)
        pos[i,0] = pos[i,0] + F1*np.cos(phi1)*dt + F2*np.cos(phi2)*dt
        pos[i,1] = pos[i,1] + F1*np.sin(phi1)*dt + F2*np.sin(phi2)*dt

    return e2e_distance(pos[0,:], pos[n-1,:])


@jit(fastmath=True, nopython=True, parallel = True)
def yes(array1, array2, time, step, length):
    L = array1
    L1 = array2
    x = np.zeros(time)
    for i in prange(L.shape[0]):
        n = 1000 + i
        pos = np.zeros((n,2))   # position array
        
        f = j = 0
        # Setting the IC
        for f in prange(n):
            pos[f,0] = f * length
        L1[i] = e2e_distance(pos[0,:], pos[n-1,:])
            
        for j in range(time):
            x[j] = animate(pos, n, length, step, 0.1)

        L[i] = np.sum(x)/time
    return L

n = 100
l = 1                  # bond length
t = int(1e10)           # simulation time
s = 1e-3               # step length
k = 500                # spring constant
L = np.zeros(n)        
L1 = np.zeros(n)       # shows initial E2E dist

start = time.time()
y = yes(L, L1, t, s, l)
x = np.linspace(10, len(y)+10, len(y))
plt.plot(x, y, 'b')
plt.plot(x, L1, 'r')
plt.xlabel('Length N')
plt.ylabel('End to End distance time average')
plt.title('Time Step' + str(s) + " Time" + str(t)+"k="+str(k))
plt.tight_layout()

l = time.localtime()
d8tym = (str(l.tm_mday) + "-" +  str(l.tm_mon) + "at" + str(l.tm_hour) + ":" + str(l.tm_min) + ":" + str(l.tm_sec))

plt.savefig("Plots/Fig_n="+str(n)+"t="+str(t)+"s="+str(s)+"k="+str(k)+"time:" + d8tym +".png")
plt.close()
end = time.time()
print("Elapsed (mins) = %3.4f" % (end/60 - start/60))
