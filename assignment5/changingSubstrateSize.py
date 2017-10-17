# File Name: changingSubstrateSize.py
# Description: See when increasing the substrate size stops giving appreciable
#               gains in accuracy vs increases in computation time
# Author: Christopher Parker
# Created: Tue Oct 10, 2017 | 01:32P EDT
# Last Modified: Tue Oct 17, 2017 | 11:15P EDT

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

import numpy as np
import scipy as sci
from numpy.linalg import norm
import time
from IPython import embed
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# start tracking the runtime of the program
start_time = time.time()

# first, we set the constants:
w = 1
sigma = 1

h_x = 1
h_y = 1

# initialize substrate size at a single atom
l = 0

# initialize vars used in computing error
norm_VDW_force = 1
final_VDW_force = 0

# this is the function that will be passed to the ODE solver as the RHS
def vdwForce(r,t):

    # compute the offsets
    k_x = (r[0]+.5)%h_x
    k_y = (r[1]+.5)%h_y

    # set the number of surrounding substrate atoms to consider
    a = np.arange(-l,l+1,1)

    # initialize total_VDW_force
    total_VDW_force = np.zeros((3,))

    # compute the distances of the floating atom from the substrate atoms
    for j in range(len(a)):
        for i in range(len(a)):

            # compute r_hat(r) = ||r||
            dx = k_x - .5 + a[i]*h_x
            dy = k_y - .5 + a[j]*h_y
            dz = r[2]
            r_hat = np.sqrt(dx**2 + dy**2 + dz**2)

            # this is the gradient of V (computed by hand). 
            gradV_common = (12*w*((sigma**6)/(r_hat**8)-(sigma**12)/(r_hat**14)))
            VDW_force = -gradV_common*np.array([dx, dy, dz])

            total_VDW_force += VDW_force

    n = len(total_VDW_force)
    final_VDW_force = total_VDW_force[n-1]

    # compute the norm of VDW Force at the final position of the particle (to
    # gauge the accuracy for each substrate size). We would like it to be
    # less than 10e-10
    global norm_VDW_force
    norm_VDW_force = norm(final_VDW_force)

    # need this statement because the first time through the while loop gives
    # a very small value of the norm (on the order of 10e-13), but then the
    # norm increases on the next step.
    if (norm_VDW_force < 10e-10):
        norm_VDW_force = 1

    return total_VDW_force

# define the time interval for odeint
t = np.linspace(0,100,501)

# define the starting point of the floater
r0 = np.array([.5, .5, .6])

# define the error tolerance for the norm of the VDW force
tol = 1.62e-6

while (norm_VDW_force >= tol):

    loop_start = time.time()
    # compute the path of the particle (and more importantly it's final resting
    # position)
    gFlow = odeint(vdwForce, r0, t, rtol=1.4e-20)
    l += 1
    m = len(gFlow)

    #print(gFlow[m-1])

    # calculate the time it takes for each run of odeint with different
    # substrate sizes
    loop_time = time.time() - loop_start

    print([norm_VDW_force, l-1, loop_time])

# stop tracking the runtime of the program
stop_time = time.time()
comp_time = stop_time - start_time

print(comp_time)
