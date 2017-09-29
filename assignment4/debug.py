# File Name: debug.py
# Description: Trying various plots for debugging gradientFlow.py
# Author: Christopher Parker
# Created: Thu Sep 28, 2017 | 12:05P EDT
# Last Modified: Fri Sep 29, 2017 | 01:37P EDT

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
from gradientFlow import vdwForce
import matplotlib as mpl

mpl.use('agg')

import matplotlib.pyplot as plt

gFlow = np.loadtxt("gradientFlow_output.txt")

# although t isn't used for anything, it must be initialized and passed
t = np.linspace(0,10,101)

# create a linspace for x to move from -.5 to .5
deltaX = np.linspace(-.5,.5,100)

# create the vector r (which represents the floating atom) w/ constant y & z
r = [deltaX,0,1]

# initialize an array to hold the values of the VDW force at each x
vdwForces = []

# now we loop through the changing x values, calling vdwForce at each step
# in order to find the total VDW force exerted on the floating atom at
# each x value
for i in range(len(deltaX)):
    vdwForces.append(vdwForce([r[0][i],r[1],r[2]],t)[0])

# now we plot deltaX on the x-axis and vdwForces on the y-axis
plt.figure()

plt.plot(deltaX,vdwForces,'k')

# label the axes
plt.xlabel('x')
plt.ylabel('Total VDW Force')

plt.savefig('deltaX_vs_VDWforce.png')
#plt.show()
