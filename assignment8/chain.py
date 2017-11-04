# File Name: chain.py
# Description: A chain of bonded atoms floating above a rectangular substrate.
# Author: Christopher Parker
# Created: Fri Nov 03, 2017 | 10:32P EDT
# Last Modified: Fri Nov 03, 2017 | 11:23P EDT

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
from numpy.linalg import norm
from scipy.integrate import odeint

if __name__ == '__main__':

    # set the constants
    w = 1
    sigma = 1
    
    h_x = 1
    h_y = 1
    
    k_s = .1
    l = 2
    
    # define the time interval for the gradient flow
    t = np.linspace(0,1,1000)
    
    # define the starting point of the floaters
    floaters = [[0 for x in range(3)] for y in range(3)]
    
    for i in range(len(floaters)):
        floaters[i][0] = i
        floaters[i][1] = i
        floaters[i][2] = 1
    
    #odeint(spring_and_VDWForce, floaters, t, atol=1.4e-10)

def spring_and_VDWForce(floaters):

    # initialize arrays for the spring forces and VDW forces
    spring_forces = [[0 for x in range(len(floaters))] for y in range(len(floaters)-1)]
    VDW_forces = np.zeros(len(floaters))

    # initialize arrays for offsets, dx, dy, dz and rhat
    k_x = np.zeros(len(floaters))
    k_y = np.zeros(len(floaters))

    dx = np.zeros(len(floaters))
    dy = np.zeros(len(floaters))
    dz = np.zeros(len(floaters))

    rhat = np.zeros(len(floaters))

    # compute offsets
    for i in range(len(k_x)):
        k_x[i] = floaters[i][0]%h_x
        k_y[i] = floaters[i][1]%h_y

    # define the size of the substrate surrounding the floaters
    a = np.arange(-2,3,1)



spring_and_VDWForce(floaters)

