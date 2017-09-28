# File Name: VDWchange_plot.py
# Description: Plot the change in VDW energy at each point along the path of
#               the floating particle
# Author: Christopher Parker
# Created: Sat Sep 23, 2017 | 11:47P EDT
# Last Modified: Wed Sep 27, 2017 | 01:51P EDT

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
import matplotlib.pyplot as plt

gFlow = np.loadtxt("gradientFlow_output.txt")


# Set the values of w and sigma
w = 1
sigma = 1

# Solve for the distances r_j
#
# First, I define the distance between the fixed atoms (assuming that there
# is one at x = 0
h_x = 1
h_y = 1


# function to compute total VDW energy using truncated sum
def VDW(x,y,z):

    # Offset k (x mod h for the x-coordinate of the floating atom)
    k_x = x%h_x
    k_y = y%h_y

    E_v = 0
    a = np.arange(-6,7,1) # set the values of r_j which we need

    for i in range(len(a)):
        for j in range(len(a)):
            r = np.sqrt((k_x + a[i]*h_x)**2 + (k_y + a[j]*h_y)**2 + z**2)
            E_v += w*((sigma**12/r**12) - 2*(sigma**6/r**6))

    return E_v

# greate an empty array to store the values of the VDW energy
VDW_energies = []
deltaVDW = [0]
z = []
y = []
x = []
# compute the total VDW energy for each point in gFlow
for ind1 in range(len(gFlow)):
    VDW_energies.append(VDW(gFlow[ind1][0], gFlow[ind1][1], gFlow[ind1][2]))

    z.append(gFlow[ind1][2])
    y.append(gFlow[ind1][1])
    x.append(gFlow[ind1][0])
    if ind1 > 0:
        # create an array of changes in the VDW energy as the floater moves
        deltaVDW.append(VDW_energies[ind1] - VDW_energies[ind1 - 1])

# plot the VDW energy of the particle over time
t = np.linspace(0,5,501)

plt.figure()
ax = plt.gca()
ax.set_xlim([0,5])
plt.plot(t,VDW_energies,'k')
plt.xlabel("Time")
plt.ylabel("VDW Energy")
plt.savefig("VDW_energies2.png")
plt.show()
