# File Name: gradientFlow.py
# Description: This code will use grad(V) to compute numerical solutions to
#               the ODEs created by the gradient flow equation.
# Author: Christopher Parker
# Created: Wed Sep 20, 2017 | 12:58P EDT
# Last Modified: Thu Oct 26, 2017 | 11:39P EDT

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
from IPython import embed
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# first, we set the constants:
w = 1
sigma = 1

h_x = 1
h_y = 1

# this is the function that will be passed to the ODE solver as the RHS
def vdwForce(r,t):

    # initialize total_VDW_force
    total_VDW_force = np.zeros((3,))

    # compute the offsets
    k_x = (r[0]+.5)%h_x
    k_y = (r[1]+.5)%h_y


    # set the number of surrounding substrate atoms to consider
    a = np.arange(-6,7,1)

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

    #print('F: ',total_VDW_force)
    return total_VDW_force


#deltaY = np.linspace(-.5,.5,10000)

# define the time interval for the gradient flow
t = np.linspace(0,10,501)

# define the starting point of the floater
r0 = np.array([.5, 2.5, .6])

# compute the gradient flow equation for each value of r, and save
# the values in an array
gFlow = odeint(vdwForce,r0,t,rtol=1.4e-20)

# write the output to a file in order to easily read and plot it

#np.savetxt("gradientFlow_output.txt",gFlow)


print(gFlow)
