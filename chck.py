import numpy as np
import scipy as sp
import random
from itertools import count
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from numba import prange,jit

box = 10
def multiple(n):
    array = np.zeros([n,2])
    plt.clf()
    plt.axis(( -1*box, box, -1*box, box))
    for i in range(len(array)):
        array[i,0] += random.randint(0,1)
        array[i,1] += random.randint(0,1)
        plt.plot(array[i,0], array[i,1], 'ro', markersize=2)
def animate(i):
    multiple(1)

ani = FuncAnimation(plt.gcf(), animate, interval=100)
plt.tight_layout()
plt.show()
