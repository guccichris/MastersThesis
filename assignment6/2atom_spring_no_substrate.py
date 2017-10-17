# File Name: 2atom_spring.py
# Description: This script will model two bonded atoms floating in space by
#               connecting them with a spring
# Author: Christopher Parker
# Created: Mon Oct 16, 2017 | 12:35P EDT
# Last Modified: Tue Oct 17, 2017 | 12:08P EDT

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
import matplotlib.pyplot as plt

# first, we set the constants:
w = 1
sigma = 1

h_x = 1
h_y = 1

k_s = 1
l = 2

# this is the function that will be passed to the ODE solver as the RHS
def springForce(r,t):

    r1 = r[:3]
    r2 = r[3:]
    # initialize total_spring_force
    total_spring_force = np.zeros((6,))

    # compute the offsets
    #k_x = (r[0]+.5)%h_x
    #k_y = (r[1]+.5)%h_y


    # set the number of surrounding substrate atoms to consider
    #a = np.arange(-2,3,1)

    # compute the distances of the floating atom from the substrate atoms
    #for j in range(len(a)):
        #for i in range(len(a)):

            # compute r_hat(r) = ||r||
            #dx = k_x - .5 + a[i]*h_x
            #dy = k_y - .5 + a[j]*h_y
            #dz = r[2]
            #r_hat = np.sqrt(dx**2 + dy**2 + dz**2)

    # this is the gradient of E_s (computed by hand). 
    gradE_r1_common = k_s*(norm(r1-r2) - l)/norm(r1 - r2)
    r1_spring_force = -gradE_r1_common*np.array([(r1[0]-r2[0]),(r1[1]-r2[1]),(r1[2]-r2[2])])
    r2_spring_force = -1*r1_spring_force

    total_spring_force = np.concatenate([r1_spring_force, r2_spring_force])

    return total_spring_force


# define the time interval for the gradient flow
t = np.linspace(0,1000,501)

# define the starting point of the floaters
r1 = np.array([1, 1, 1])
r2 = np.array([0, 0, 1])
r0 = np.concatenate([r1, r2])

# compute the gradient flow equation for each value of t, and save
# the values in an array
gFlow = odeint(springForce,r0,t,rtol=1.4e-20)

for elt in gFlow:
    r1 = elt[:3]
    r2 = elt[3:]
    dist = norm(r1 - r2)
    print(r1, r2)
    print(dist)

print(gFlow)
