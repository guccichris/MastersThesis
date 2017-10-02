# File Name: gradientFlow.py
# Description: This code will use grad(V) to compute numerical solutions to
#               the ODEs created by the gradient flow equation.
# Author: Christopher Parker
# Created: Wed Sep 20, 2017 | 12:58P EDT
# Last Modified: Mon Oct 02, 2017 | 02:07P EDT

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
    k_x = r[0]%h_x
    k_y = r[1]%h_y


    # set the number of surrounding substrate atoms to consider
    a = np.arange(-6,7,1)

    # compute the distances of the floating atom from the substrate atoms
    for j in range(len(a)):
        for i in range(len(a)):

            # compute r_hat(r) = ||r||
            dx = k_x + a[i]*h_x
            dy = k_y + a[j]*h_y
            dz = r[2]
            r_hat = np.sqrt(dx**2 + dy**2 + dz**2)

            #r2 = [a[i],a[j],0]
            #r_hat = dist(r,r2)
            #d = r - r2
            # this is the gradient of V (computed by hand). 
            gradV_common = (12*w*((sigma**6)/(r_hat**8)-(sigma**12)/(r_hat**14)))
            VDW_force = -gradV_common*np.array([dx, dy, dz])

            total_VDW_force += VDW_force


    #print('total_VDW_force',total_VDW_force)
    return total_VDW_force

def dist(r1,r2):
    return np.sqrt((r1[0]-r2[0])**2 + (r1[1] - r2[1])**2 + (r1[2] - r2[2])**2)

# define the time interval for the gradient flow
t = np.linspace(0,100,5001)

# define the starting point of the floater
r0 = np.array([0, 0, .6])
print(vdwForce(r0,t))

# compute the gradient flow equation for each value of r, and save
# the values in an array
gFlow = odeint(vdwForce,r0,t,rtol=1.4e-14)

# write the output to a file in order to easily read and plot it

#np.savetxt("gradientFlow_output.txt",gFlow)


print(gFlow)
