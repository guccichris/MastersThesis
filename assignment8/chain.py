# File Name: chain.py
# Description: A chain of bonded atoms floating above a rectangular substrate.
# Author: Christopher Parker
# Created: Fri Nov 03, 2017 | 10:32P EDT
# Last Modified: Mon Nov 06, 2017 | 01:53P EST

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
from IPython import embed
import pprint

# define the number of floaters
floaters = np.zeros(4*3)
n = len(floaters)

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
    gradV_common = np.zeros(m)
    VDW_force = [[0 for i in range(3)] for j in range(m)]
    spring_force = [[0 for i in range(3)] for j in range(m)]

    # initialize arrays for offsets, dx, dy, dz and rhat
    k_x = np.zeros(m)
    k_y = np.zeros(m)

    dx = np.zeros(m)
    dy = np.zeros(m)
    dz = np.zeros(m)

    # define the size of the substrate surrounding the floaters
    a = np.arange(-6,7,1)

    # rhat needs to be a 3x5x5 matrix in this case
    rhat = np.zeros(m)

    # compute offsets
    for i in range(m):
        k_x[i] = (floaters[i*3] + .5)%h_x
        k_y[i] = (floaters[i*3] + .5)%h_y

    # compute the distances of the floating atom from the substrate atoms
    # and the VDW forces acting on them as a result
    for j in range(len(a)):
        for i in range(len(a)):
            for k in range(m):

                # compute dx, dy, dz and rhat
                dx[k] = k_x[k] - .5 + a[i]*h_x
                dy[k] = k_y[k] - .5 + a[j]*h_y
                dz[k] = floaters[k*3 + 2]
                rhat[k] = np.sqrt(dx[k]**2 + dy[k]**2 + dz[k]**2)

                # this is the gradient of V (computed by hand), it is used to compute the
                # VDW forces acting on the floating atoms
                gradV_common[k] = (12*w*((sigma**6)/(rhat[k]**8)-(sigma**12)/(rhat[k]**14)))
                VDW_force[k] = -gradV_common[k]*np.array([dx[k], dy[k], dz[k]])

            # combine the VDW forces for all floating atoms
            combined_VDW_forces += np.hstack(VDW_force)

    # this is the gradient of E_s (computed by hand), it is used to compute the
    # spring forces acting on the atoms
    for i in range(m-1):
        gradE_common = k_s*(norm(r[i]-r[i+1]) - l)/norm(r[i] - r[i+1])
        spring_force[i] += -gradE_common*np.array([(r[i][0]-r[i+1][0]),(r[i][1]-r[i+1][1]),(r[i][2]-r[i+1][2])])
        spring_force[i+1] = gradE_common*np.array([(r[i][0]-r[i+1][0]),(r[i][1]-r[i+1][1]),(r[i][2]-r[i+1][2])])


    combined_spring_forces = np.hstack(spring_force)
    total_spring_and_VDW_force = np.add(combined_spring_forces, combined_VDW_forces)

    return total_spring_and_VDW_force

if __name__ == '__main__':

    # set the constants
    w = 1
    sigma = 1

    h_x = 1
    h_y = 1

    k_s = 1
    l = 1

    # define the time interval for the gradient flow
    t = np.linspace(0,1,1000)

    # set the positions of the floating atoms
    count = 0
    for i in range(0,n,3):
        floaters[i] = count
        floaters[i+1] = count
        floaters[i+2] = 1
        count += 1

    # compute the changing positions of the atoms as a result of the spring
    # and VDW forces acting on them
    gFlow = odeint(spring_and_VDWForce, floaters, t, atol=1.4e-10)

    # write the results to gFlow_chain.txt for use in plotting
    np.savetxt('gFlow_chain.txt', gFlow)

    pprint.pprint(gFlow)
    pprint.pprint(gFlow[-1])
