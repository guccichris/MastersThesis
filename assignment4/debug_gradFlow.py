# File Name: debug_gradFlow.py
# Description: Trying the looping in gradientFlow.py with a known function
#               to ensure that I am getting the proper output
# Author: Christopher Parker
# Created: Fri Sep 29, 2017 | 01:54P EDT
# Last Modified: Tue Oct 03, 2017 | 11:06P EDT

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
from IPython import embed
from scipy.integrate import odeint

# I am just going to copy over the code from gradientFlow.py except with a
# known (simple) system of ODEs
#
# (Just have to keep in mind that I need to keep the nested for loops in order
# to bebug properly)

h_x = 1
h_y = 1

# this is the function that will be passed to the ODE solver as the RHS
def simpleODE(r,t):

    # compute the offsets
    k_x = r[0]%h_x
    k_y = r[1]%h_y

    # need to loop through a 2x2 matrix for w/e this function is going to be
    #
    # set the number of surrounding substrate atoms to consider
    a = np.arange(-6,7,1)

    # compute the distances of the floating atom from the substrate atoms
    for j in range(len(a)):
        for i in range(len(a)):

            # compute r_hat(r) = ||r||
            dx = k_x + a[i]*h_x
            dy = k_y + a[j]*h_y
            dz = r[2]
            print('dx: ', dx, 'dy: ', dy, 'dz: ', dz)
            r_hat = np.sqrt(dx**2 + dy**2 + dz**2)
            print('r_hat: ', r_hat)


    #return RETURN_VALUE  # not sure what this will be yet


# define the time interval for the gradient flow
t = np.linspace(0,100,5001)

# define the starting point of the floater
r0 = np.array([1, 1, 1])

# compute the gradient flow equation for each value of r, and save
# the values in an array
gFlow = odeint(simpleODE,r0,t,rtol=1.4e-20)


print(gFlow)
