# File Name: debug.py
# Description: Trying various plots for debugging gradientFlow.py
# Author: Christopher Parker
# Created: Thu Sep 28, 2017 | 12:05P EDT
# Last Modified: Mon Oct 30, 2017 | 07:08P EDT

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
#import matplotlib as mpl

#mpl.use('agg')

import matplotlib.pyplot as plt

# although t isn't used for anything, it must be initialized and passed
t = np.linspace(0,10,101)

# create a linspace for x to move from -.5 to .5 (and same for y)
deltaX = np.linspace(-.5,.5,10000)
deltaY = np.linspace(-.5,.5,10000)
# create a linspace for z to move from .8 to 1.2
deltaZ = np.linspace(.8,1.2,1000)

# create the vector r (which represents the floating atom) w/ constant y & z
r_x = [deltaX, 0, 1]

# this time create r with z changing, x & y constant
r_z = [0, 0, deltaZ]

# this time create r with y changing, x & z constant
r_y = [1, deltaY, .8]

# initialize an array to hold the values of the VDW force at each x
vdwForces = []

# now we loop through the changing x values, calling vdwForce at each step
# in order to find the total VDW force exerted on the floating atom at
# each x value
for i in range(len(deltaZ)):
    vdwForces.append(vdwForce([r_z[0],r_z[1],r_z[2][i]],t)[2])

# now we plot deltaX on the x-axis and vdwForces on the y-axis
plt.figure()

plt.plot(deltaZ,vdwForces,'k')

plt.plot([.8, 1.2],[0, 0])

# label the axes
plt.xlabel('z')
plt.ylabel('Total VDW Force')

# set axes limits
#plt.xlim([-1,1])
#plt.ylim([-1,1])

#plt.savefig('deltaY_vs_VDWforce_5_8.png')
plt.show()
