# File Name: 2atom_spring.py
# Description: This script will model two bonded atoms floating in space by
#               connecting them with a spring
# Author: Christopher Parker
# Created: Mon Oct 16, 2017 | 12:35P EDT
# Last Modified: Sat Nov 04, 2017 | 02:06P EDT

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
from IPython import embed
from scipy.integrate import odeint


# first, we set the constants:
w = 1
sigma = 1

h_x = 1
h_y = 1

k_s = .1
l = 2

# this is the function that will be passed to the ODE solver as the RHS
def spring_and_VDWForce(r,t):

    # split r1 and r2
    r1 = r[:3]
    r2 = r[3:]

    # initialize total_spring_force and total_VDW_force
    total_spring_force = np.zeros(6)
    total_VDW_force = np.zeros(6)

    # initialize offsets, dx, dy, dz and r_hat as arrays
    k_x = np.zeros(2)
    k_y = np.zeros(2)
    dx = np.zeros(2)
    dy = np.zeros(2)
    dz = np.zeros(2)
    r_hat = np.zeros(2)

    # compute the offsets
    k_x[0] = (r1[0]+.5)%h_x
    k_y[0] = (r1[1]+.5)%h_y

    k_x[1] = (r2[0]+.5)%h_x
    k_y[1] = (r2[1]+.5)%h_y

    # set the number of surrounding substrate atoms to consider
    a = np.arange(-2,3,1)

    # compute the distances of the floating atom from the substrate atoms
    for j in range(len(a)):
        for i in range(len(a)):

            # compute r_hat(r) = ||r||
            dx[0] = k_x[0] - .5 + a[i]*h_x
            dy[0] = k_y[0] - .5 + a[j]*h_y
            dz[0] = r1[2]
            r_hat[0] = np.sqrt(dx[0]**2 + dy[0]**2 + dz[0]**2)

            dx[1] = k_x[1] - .5 + a[i]*h_x
            dy[1] = k_y[1] - .5 + a[j]*h_y
            dz[1] = r2[2]
            r_hat[1] = np.sqrt(dx[1]**2 + dy[1]**2 + dz[1]**2)

            # this is the gradient of V (computed by hand). 
            gradV_common_r1 = (12*w*((sigma**6)/(r_hat[0]**8)-(sigma**12)/(r_hat[0]**14)))
            VDW_force_r1 = -gradV_common_r1*np.array([dx[0], dy[0], dz[0]])

            gradV_common_r2 = (12*w*((sigma**6)/(r_hat[1]**8)-(sigma**12)/(r_hat[1]**14)))
            VDW_force_r2 = -gradV_common_r2*np.array([dx[1], dy[1], dz[1]])

            total_VDW_force += np.concatenate([VDW_force_r1, VDW_force_r2])

    # this is the gradient of E_s (computed by hand). 
    gradE_r1_common = k_s*(norm(r1-r2) - l)/norm(r1 - r2)
    r1_spring_force = -gradE_r1_common*np.array([(r1[0]-r2[0]),(r1[1]-r2[1]),(r1[2]-r2[2])])
    r2_spring_force = -1*r1_spring_force

    total_spring_force = np.concatenate([r1_spring_force, r2_spring_force])
    total_spring_and_VDW_force = np.add(total_spring_force, total_VDW_force)

    #print(total_VDW_force, total_spring_force)

    return total_spring_and_VDW_force


# define the time interval for the gradient flow
t = np.linspace(0,1,1000)

# define the starting point of the floaters
r1 = np.array([1, 1, 1])
r2 = np.array([0, 0, 1])
r0 = np.concatenate([r1, r2])

# compute the gradient flow equation for each value of t, and save
# the values in an array
gFlow = odeint(spring_and_VDWForce,r0,t,rtol=1.4e-1)

print(gFlow)
np.savetxt('gFlow_2atoms.txt', gFlow)
