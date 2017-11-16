# File Name: total_energies.py
# Description: Compute the total energies acting on 2 bonded atoms over a
#               substrate
# Author: Christopher Parker
# Created: Tue Oct 31, 2017 | 11:32P EDT
# Last Modified: Mon Nov 13, 2017 | 07:59P EST

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
import matplotlib.pyplot as plt

# first, we set the constants:
w = 1
sigma = 1

h_x = 1
h_y = 1

k_s = 0
l = 2

from chain import VDW, spring_and_VDWForce
def spring(r1,r2):
    E_s = .5*k_s*(norm(r1 - r2) - l)**2
    return E_s

# define the time interval for the gradient flow
t = np.linspace(0,1,1000)

# define the starting point of the floaters
r1 = np.array([1, 1, 1])
r2 = np.array([0, 0, 1])
r0 = np.concatenate([r1, r2])

# compute the gradient flow equation for each value of t, and save
# the values in an array
gFlow = odeint(spring_and_VDWForce,r0,t,rtol=1.4e-1)

floater1 = gFlow[:,:3]
floater2 = gFlow[:,3:]

VDW_energy1 = 0
VDW_energy2 = 0
spring_energy = 0
total_energy = np.zeros(len(floater1))
TE_floater1 = np.zeros(len(floater1))
TE_floater2 = np.zeros(len(floater2))
delta_TE = np.zeros(len(floater1)-1)

for i in range(len(floater1)):
    VDW_energy1 = VDW(floater1[i])
    VDW_energy2 = VDW(floater2[i])
    spring_energy = spring(floater1[i],floater2[i])

    total_energy[i] = VDW_energy1 + VDW_energy2 + spring_energy
    TE_floater1[i] = VDW_energy1 + spring_energy
    TE_floater2[i] = VDW_energy2 + spring_energy
    if i > 0:
        delta_TE[i-1] = total_energy[i] - total_energy[i-1]

#for i in range(len(delta_TE)):
    #if delta_TE[i] >= 0:
        #print(i+1)
        #print(delta_TE[i])
        #print(gFlow[i+1])

#count = 0
#for i in range(len(TE_floater1)):
    #if TE_floater1[i] == TE_floater2[i]:
        #count += 1
    #else:
        #print(TE_floater1[i],TE_floater2[i],sep=', ')

#print(count, len(TE_floater1))

# create a figure with 4 subplots
f, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True)

# plot the total energy of the system and each individual atom
ax1.plot(t, total_energy, 'k', label="Total")
ax1.plot(t, TE_floater1, 'b', label="Atom 1")
ax1.plot(t, TE_floater2, 'r', label="Atom 2")
ax1.set_title("Energies")
ax1.legend()

# plot the x coords of both atoms
ax2.plot(t, gFlow[:,0], 'b', label="Atom 1")
ax2.plot(t, gFlow[:,3], 'r', label="Atom 2")
ax2.set_title("x coordinates")
ax2.legend()

# plot the y coords of both atoms
ax3.plot(t, gFlow[:,1], 'b', label="Atom 1")
ax3.plot(t, gFlow[:,4], 'r', label="Atom 2")
ax3.set_title("y coordinates")
ax3.legend()

# plot the z coords of both atoms
ax4.plot(t, gFlow[:,2], 'b', label="Atom 1")
ax4.plot(t, gFlow[:,5], 'r', label="Atom 2")
ax4.set_title("z coordinates")
ax4.legend()

# create individual plot for total energy
f2 = plt.figure()
ax = f2.add_subplot(111)
ax.plot(t, total_energy, 'k')

plt.show()
