# File Name: chain.py
# Description: A chain of bonded atoms floating above a rectangular substrate.
# Author: Christopher Parker
# Created: Fri Nov 03, 2017 | 10:32P EDT
# Last Modified: Wed Nov 22, 2017 | 01:51P EST

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#                           GNU GPL LICENSE                            #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
#                                                                      #
# Copyright Christopher Parker 2017 <cjp65@case.edu>                   #
#                                                                      #
# This program is free software: you can redistribute it and/or modify #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation, either version 3 of the License, or    #
# (at your option) any later version.                                  #
#                                                                      #
# This program is distributed in the hope that it will be useful,      #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         #
# GNU General Public License for more details.                         #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with this program. If not, see <http://www.gnu.org/licenses/>. #
#                                                                      #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm
from scipy.integrate import odeint
from IPython import embed
import pprint
import time

start_time = time.time()
# define the number of floaters
floaters = np.zeros(10*3)
n = len(floaters)
m = n/3

# set number of substrate atoms (substrate will be 2*ns + 1 side lengths)
ns = 5

# this is the RHS of the system of ODEs solved by odeint. It returns the forces
# acting on each floater in the x, y and z directions
def spring_and_VDWForce(floaters, t):

    # compute number of atoms
    m = int(len(floaters)/3)
    n = int(len(floaters))

    # split the atoms into 1x3 arrays
    r = [[0 for i in range(3)] for j in range(m)]
    for i in range(m):
        r[i] = floaters[i*3:i*3+3]

    # initialize arrays for the spring forces and VDW forces
    combined_spring_forces = np.zeros(n)
    combined_VDW_forces = np.zeros(n)
    VDW_force = [[0 for i in range(3)] for j in range(m)]
    spring_force = [[0 for i in range(3)] for j in range(m)]

    # define the size of the substrate surrounding the floaters
    a = np.arange(-ns,ns+1,1)

    # compute the distances of the floating atom from the substrate atoms
    # and the VDW forces acting on them as a result
    for j in range(len(a)):
        for i in range(len(a)):
            for k in range(m):

                # compute the offsets
                k_x = (r[k][0] + .5)%h_x
                k_y = (r[k][1] + .5)%h_y

                # compute dx, dy, dz and rhat
                dx = k_x - .5 + a[i]*h_x
                dy = k_y - .5 + a[j]*h_y
                dz = r[k][2]
                rhat = np.sqrt(dx**2 + dy**2 + dz**2)

                # this is the gradient of V (computed by hand), it is used to compute the
                # VDW forces acting on the floating atoms
                gradV_common = (12*w*((sigma**6)/(rhat**8)-(sigma**12)/(rhat**14)))
                VDW_force[k] = -gradV_common*np.array([dx, dy, dz])

            # combine the VDW forces for all floating atoms
            combined_VDW_forces += np.hstack(VDW_force)

    # this is the gradient of E_s (computed by hand), it is used to compute the
    # spring forces acting on the atoms
    for i in range(m-1):
        gradE_common = k_s*(norm(r[i]-r[i+1]) - l)/norm(r[i] - r[i+1])
        current_sf = -gradE_common*np.array([(r[i][0]-r[i+1][0]),(r[i][1]-r[i+1][1]),(r[i][2]-r[i+1][2])])
        spring_force[i] += current_sf
        spring_force[i+1] = -current_sf


    combined_spring_forces = np.hstack(spring_force)
    total_spring_and_VDW_force = np.add(combined_spring_forces, combined_VDW_forces)

    return total_spring_and_VDW_force

# function to compute VDW energy for plotting
def VDW(r):

    x = r[0]
    y = r[1]
    z = r[2]

    # Offset k (x mod h for the x-coordinate of the floating atom)
    k_x = x%h_x
    k_y = y%h_y

    E_v = 0
    a = np.arange(-ns,ns+1,1) # set the values of r_j which we need

    for i in range(len(a)):
        for j in range(len(a)):
            r = np.sqrt((k_x + a[i]*h_x)**2 + (k_y + a[j]*h_y)**2 + z**2)
            E_v += w*((sigma**12/r**12) - 2*(sigma**6/r**6))

    return E_v

if __name__ == '__main__':

    # set the constants
    w = 1
    sigma = 1

    h_x = 1
    h_y = 1

    k_s = 5
    l = 1

    # define the time interval for the gradient flow
    t = np.linspace(0,6,20)

    # set the positions of the floating atoms
    count = 0
    for i in range(0,n,3):
        floaters[i] = count
        floaters[i+1] = 1.7*count
        floaters[i+2] = 1
        count += 1

    # compute the changing positions of the atoms as a result of the spring
    # and VDW forces acting on them
    gFlow = odeint(spring_and_VDWForce, floaters, t, atol=1.4e-10)

    # compute and print the distance between each pair of bonded atoms (for
    # a 4 atom chain) for debugging purposes
    #dist1 = norm(gFlow[-1,:3]-gFlow[-1,3:6])
#    dist2 = norm(gFlow[-1,3:6]-gFlow[-1,6:9])
#    dist3 = norm(gFlow[-1,6:9]-gFlow[-1,9:])
#    print(dist1,dist2,dist3, sep=', ')

    # write the results to gFlow_chain.txt for use in plotting
    np.savetxt('gFlow_chain.txt', gFlow)

    # print the output of odeint
    pprint.pprint(gFlow)
    pprint.pprint(gFlow[-1])

    end_time = time.time()
    print(end_time - start_time)
