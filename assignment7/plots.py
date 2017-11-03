# File Name: plots.py
# Description: This will create plots of distance between the two atoms vs time,
#               paths of the 2 atoms vs time and position of the center of mass
#               of the 2 atoms vs time
# Author: Christopher Parker
# Created: Mon Oct 30, 2017 | 09:30P EDT
# Last Modified:

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
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from numpy.linalg import norm

# import the path data from 2atom_spring.py
gFlow = np.loadtxt('gFlow_2atoms.txt')

# separate the 2 floaters in the array gFlow
floater1 = gFlow[:,:3]
floater2 = gFlow[:,3:]

# define the time interval over which the plots range
t = np.linspace(0,1,1000)

#
# plot the distance between the atoms vs time
#

dist_between = np.zeros(len(floater1))

for i in range(len(floater1)):
    dist_between[i] = norm(floater1[i] - floater2[i])

# create the plot
fig = plt.figure(1)
ax = fig.add_subplot(111)
ax.plot(t, dist_between, 'k')

# set the axis labels
ax.set_xlabel('Time')
ax.set_ylabel('Distance Between Floating Atoms')

# save the plot as dist_vs_time.png
plt.savefig('dist_vs_time.png')

#
# plot the paths of the atoms over time
#

# create 3d figure
fig = plt.figure()
ax = fig.gca(projection='3d')

# plot position of floater1 vs time
ax.plot(floater1[:,0],floater1[:,1],floater1[:,2], label='floater1')

# plot position of floater2 vs time
ax.plot(floater2[:,0], floater2[:,1], floater2[:,2], label='floater2')

# show legend
ax.legend()

# save plot to paths_of_atoms.png
plt.savefig('paths_of_atoms.png')

#
# plot position of center of mass of floating atoms vs time
#

# since we are assuming the floating atoms have equal mass, the center of mass
# of 2 floating atoms is simply the midpoint between them
midpt = np.zeros(len(floater1))

midpt = (floater1 + floater2)/2
print(midpt)

# plot the midpoint vs time
fig = plt.figure()
ax1 = fig.gca(projection='3d')
ax1.plot(midpt[:,0], midpt[:,1], midpt[:,2])

# save the plot as cm_vs_time.png
plt.savefig('cm_vs_time.png')


plt.show()
