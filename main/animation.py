import numpy as np
import random
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from numba import prange,jit

box = 10    #area of graph

pos = np.zeros([2,2])   #position array
phi = np.zeros([1])     #angle array
length = 2              #distance between the two particles
mks = 5                 #size of the two particles
pos[1,0] = 4            #arbitrary initial position
dt = .5                #delta_t
sqrt_dt = np.sqrt(dt)

@jit(fastmath=True)
def angle(position):
    x = position[1,0] - position[0,0]
    y = position[1,1] - position[0,1]
    return np.arctan2(y,x)

@jit(fastmath=True)
def spring(position):
    k = .2
    x_1 = position[0,0]
    y_1 = position[0,1]
    x_2 = position[1,0]
    y_2 = position[1,1]
    distance = np.sqrt((x_1 - x_2)**2 + (y_1 - y_2)**2)
    return -k*(distance - length)


@jit(fastmath=True)
def pos_update(position, angle):
    F = spring(position)
    for i in range(len(pos)):
        sign = -1 if (i % 2)==0 else +1
        position[i,0] += np.random.normal(0, 0.5)*sqrt_dt + sign*F*np.cos(angle)*dt
        position[i,1] += np.random.normal(0, 0.5)*sqrt_dt + sign*F*np.sin(angle)*dt
        



    #position[0,0] += (-1)**random.randint(0,1)*np.random.normal()*np.sqrt(dt) - F*np.cos(angle)*dt
    #position[0,1] += (-1)**random.randint(0,1)*np.random.normal()*0.4 - F*np.sin(angle)
    #position[1,0] += (-1)**random.randint(0,1)*np.random.normal()*0.4 + F*np.cos(angle) #+ length*np.cos(angle[0])
    #position[1,1] += (-1)**random.randint(0,1)*np.random.normal()*0.4 + F*np.sin(angle) #+ length*np.sin(angle[0])


def animate(i):
    phi = angle(pos)
    pos_update(pos, phi)
    plt.clf()
    plt.axis(( -1*box, box, -1*box, box))
    plt.plot(pos[0,0], pos[0,1], 'ro', markersize=mks)
    plt.plot(pos[1,0], pos[1,1], 'bo', markersize=mks)
    plt.plot(pos[:,0], pos[:,1], 'gray', linestyle=':', marker='')

ani = FuncAnimation(plt.gcf(), animate, interval=100)
plt.tight_layout()
plt.show()
