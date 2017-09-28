# Code to graph Van der Waals energy based on distance between atoms

import matplotlib.pyplot as plt
import numpy as np
from LJ12_6 import LJ

# Create an array for the varying distance between atoms
r = np.arange(0.01, 10., 0.01)

# For each value of r, plug it into LJ and create an array V of these values
V = []
for i in range(len(r)):
    V.append(LJ(r[i]))


# Create a graph with r on the x-axis and V on the y-axis
plt.plot(r,V,'k')

# Label the axes
plt.xlabel('r', size='x-large')
plt.ylabel('VDW energy')

# Set where grid lines will be on the graph
x_ticks = np.arange(0,11,1)
y_ticks = np.arange(-2,11,1)
plt.xticks(x_ticks)
plt.yticks(y_ticks)

# Define the range of the x-axis and y-axis, then overlay the grid
plt.axis([0,10,-2, 10])
plt.grid(True)

# Place axes at x = 0 and  y = 0 rather than at bottom and left
plt.gca().spines['left'].set_position('zero')
plt.gca().spines['bottom'].set_position('zero')

plt.savefig('VDWenergy_vs_distance.png')
plt.show()
