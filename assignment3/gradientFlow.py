# File Name: gradientFlow.py
# Description: This code will use grad(V) to compute numerical solutions to
#               the ODEs created by the gradient flow equation.
# Author: Christopher Parker
# Created: Wed Sep 20, 2017 | 12:58P EDT
# Last Modified: Sat Sep 23, 2017 | 12:07P EDT

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
from scipy.integrate import odeint

# first, we set the constants:
w = 1
sigma = 1

h_x = 1
h_y = 1

# this is the function that will be passed to the ODE solver as the RHS
def vdwForce(r,t):

    # initialize total_VDW_force
    total_VDW_force = [0,0,0]

    # compute the offsets
    k_x = r[0]%h_x
    k_y = r[1]%h_y


    # set the number of surrounding substrate atoms to consider
    a = np.arange(-6,6,1)

    # compute the distances of the floating atom from the substrate atoms
    for i in range(len(a)):
        for j in range(len(a)):

            # compute r_hat(r) = ||r||

            dx = k_x + a[i]*h_x
            dy = k_y + a[j]*h_y
            dz = r[2]
            r_hat = np.sqrt(dx**2 + dy**2 + dz**2)

            # this is the gradient of V (computed by hand). 
            gradV_common = (12*w*(-(sigma**12)/(r_hat**13)+(sigma**6)/(r_hat**7)))*(1/r_hat)
            gradV = [gradV_common*dx, gradV_common*dy, gradV_common*dz]

            # convert the x, y and z components to floats and multiply by -1
            drdx = float(-gradV[0])
            drdy = float(-gradV[1])
            drdz = float(-gradV[2])

            VDW_force = [drdx, drdy, drdz]

            total_VDW_force[0] += VDW_force[0]
            total_VDW_force[1] += VDW_force[1]
            total_VDW_force[2] += VDW_force[2]

    return total_VDW_force


# create the function which will compute the gradient flow equation given a
# floating atom
def gradFlow(x,y,z):

    # let's try it with a constant r to start
    r0 = [x, y, z]


    # compute the gradient flow equation for each value of r, and save
    # the values in an array
    gradientFlow = odeint(vdwForce,r0,t)
    return gradientFlow


# define the time interval for the gradient flow
t = np.linspace(0,10,101)

# call the function gradFlow for some floater in order to compute the path of
# the particle over time.
gFlow = gradFlow(1.5,1.9,1)

# write the output to a file in order to easily read and plot it
filename = 'gradientFlow_output.dat'
outFile = open(filename, 'w')
gFlow.tofile(outFile, sep="\n")
outFile.close()

readFile = open(filename, 'r')
test = np.fromfile(readFile)
print(test)
readFile.close()

#print(gFlow)
