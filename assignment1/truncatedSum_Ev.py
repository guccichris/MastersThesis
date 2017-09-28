# Compute the truncated sum E_v for total VDW energy
#
# now, the substrate will be a plane rather than a line

import numpy as np

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
floater = (np.linspace(1,10,100),np.linspace(1,10,100),1)
for pt in floater[0]:
    



#if type(floater[0]).__module__ == np.__name__:
#    for pt in floater[0]:
#        print(VDW(pt,floater[1],floater[2]))
#elif type(floater[1]).__module__ == np.__name__:
#    for pt in floater[1]:
#        print(VDW(floater[0],pt,floater[2]))
#elif type(floater[2]).__module__ == np.__name__:
#    for pt in floater[2]:
#        print(VDW(floater[0],floater[1],pt))
#else:
#    print(VDW(floater[0],floater[1],floater[2]))
