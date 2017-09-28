# Compute the truncated sum E_v for total VDW energy

import numpy as np

# Set the values of w and sigma
w = 1
sigma = 1

# Solve for the distances r_j
#
# First, I define the distance between the fixed atoms (assuming that there
# is one at x = 0
h = 1.5


# function to compute total VDW energy using truncated sum
def VDW(x,y):

    # Offset k (x mod h for the x-coordinate of the floating atom)
    k = x%h

    E_v = 0
    j = np.arange(-6,6,1) # set the values of r_j which we need

    for i in range(len(j)):
        r = np.sqrt((k+j[i]*h)**2 + y**2)
        E_v += w*((sigma**12/r**12) - 2*(sigma**6/r**6))

    return E_v


# create a function to compute VDW energy along a line, then I can make a script to plot this interaction
# so I created a function, now I can call it at any point (x,y). so if I create a line,
# I can call it at every point along that line.
floater = (np.linspace(1,10,100),1)
for pt in floater[0]:
    print(VDW(pt,floater[1]))
