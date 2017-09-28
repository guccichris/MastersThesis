# Compute the truncated sum E_v for total VDW energy
#
# now, the substrate will be a plane rather than a line

import numpy as np
import matplotlib as mpl

# allows matplotlib to save files in .png format
#mpl.use('agg')

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

# Set the values of w and sigma
w = 1
sigma = 1

# Solve for the distances r_j
#
# First, I define the distance between the fixed atoms (assuming that there
# is one at x = 0
h_x = 1
h_y = 1.5


# function to compute total VDW energy using truncated sum
def VDW(x,y,z):

    # Offset k (x mod h for the x-coordinate of the floating atom)
    k_x = x%h_x
    k_y = y%h_y

    E_v = 0
    a = np.arange(-6,6,1) # set the values of r_j which we need

    for i in range(len(a)):
        for j in range(len(a)):
            r = np.sqrt((k_x + a[i]*h_x)**2 + (k_y + a[j]*h_y)**2 + z**2)
            E_v += w*((sigma**12/r**12) - 2*(sigma**6/r**6))

    return E_v


# so I created a function, now I can call it at any point (x,y,z). so if I create a line,
# I can call it at every point along that line.
Z = np.zeros((100,100))
floater = (np.linspace(1,10,100),np.linspace(1,10,100),0.8)
for k in range(100):
    for j in range(100):
       Z[j,k] = VDW(floater[0][k],floater[1][j],floater[2])

# so the matrix Z are the z values we want to plot, I simply need to create
# the appropriate arrays for x and y values in order to create the 3D surface
# plot with matplotlib.

X = np.linspace(1,10,100)
Y = np.linspace(1,10,100)

# here starts the process of creating the 3d plot
fig = plt.figure()
ax = fig.gca(projection='3d')

# create a numpy meshgrid using the X and Y arrays
X, Y = np.meshgrid(X, Y)

# plot the surface
surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=True)

# add a color bar which maps values to colors
fig.colorbar(surf, shrink=0.5, aspect=5)

# set z-axis bounds
#ax.set_zlim3d(-3,0)

# label axes
ax.set_zlabel('VDW Energy')
ax.set_xlabel('X')
ax.set_ylabel('Y')

plt.show()
