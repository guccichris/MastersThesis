# File Name: total_energies.py
# Description: Compute the total energies acting on 2 bonded atoms over a
#               substrate
# Author: Christopher Parker
# Created: Tue Oct 31, 2017 | 11:32P EDT
# Last Modified: Tue Nov 07, 2017 | 01:41P EST

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

# this is the function that will be passed to the ODE solver as the RHS
def spring_and_VDWForce(r,t):

    # split r1 and r2
    r1 = r[:3]
    r2 = r[3:]

    # initialize total_spring_force and total_VDW_force
    total_spring_force = np.zeros(6)
    total_VDW_force = np.zeros(6)

    # initialize offsets, dx, dy, dz and r_hat as arrays
    k_x = np.zeros(2)
    k_y = np.zeros(2)
    dx = np.zeros(2)
    dy = np.zeros(2)
    dz = np.zeros(2)
    r_hat = np.zeros(2)

    # compute the offsets
    k_x[0] = (r1[0]+.5)%h_x
    k_y[0] = (r1[1]+.5)%h_y

    k_x[1] = (r2[0]+.5)%h_x
    k_y[1] = (r2[1]+.5)%h_y

    # set the number of surrounding substrate atoms to consider
    a = np.arange(-1,2,1)

    # compute the distances of the floating atom from the substrate atoms
    for j in range(len(a)):
        for i in range(len(a)):

            # compute r_hat(r) = ||r||
            dx[0] = k_x[0] - .5 + a[i]*h_x
            dy[0] = k_y[0] - .5 + a[j]*h_y
            dz[0] = r1[2]
            r_hat[0] = np.sqrt(dx[0]**2 + dy[0]**2 + dz[0]**2)

            dx[1] = k_x[1] - .5 + a[i]*h_x
            dy[1] = k_y[1] - .5 + a[j]*h_y
            dz[1] = r2[2]
            r_hat[1] = np.sqrt(dx[1]**2 + dy[1]**2 + dz[1]**2)


            # this is the gradient of V (computed by hand), it is used to compute the
            # VDW forces acting on the floating atoms
            gradV_common_r1 = (12*w*((sigma**6)/(r_hat[0]**8)-(sigma**12)/(r_hat[0]**14)))
            VDW_force_r1 = -gradV_common_r1*np.array([dx[0], dy[0], dz[0]])

            gradV_common_r2 = (12*w*((sigma**6)/(r_hat[1]**8)-(sigma**12)/(r_hat[1]**14)))
            VDW_force_r2 = -gradV_common_r2*np.array([dx[1], dy[1], dz[1]])

            total_VDW_force += np.concatenate([VDW_force_r1, VDW_force_r2])


    # this is the gradient of E_s (computed by hand), it is used to compute the
    # spring forces acting on the atoms
    gradE_r1_common = k_s*(norm(r1-r2) - l)/norm(r1 - r2)
    r1_spring_force = -gradE_r1_common*np.array([(r1[0]-r2[0]),(r1[1]-r2[1]),(r1[2]-r2[2])])
    r2_spring_force = -1*r1_spring_force

    total_spring_force = np.concatenate([r1_spring_force, r2_spring_force])
    total_spring_and_VDW_force = np.add(total_spring_force, total_VDW_force)

    #print(total_VDW_force, total_spring_force)

    return total_spring_and_VDW_force

# function to compute total VDW energy using truncated sum
def VDW(r):

    x = r[0]
    y = r[1]
    z = r[2]

    # Offset k (x mod h for the x-coordinate of the floating atom)
    k_x = x%h_x
    k_y = y%h_y

    E_v = 0
    a = np.arange(-1,2,1) # set the values of r_j which we need

    for i in range(len(a)):
        for j in range(len(a)):
            r = np.sqrt((k_x + a[i]*h_x)**2 + (k_y + a[j]*h_y)**2 + z**2)
            E_v += w*((sigma**12/r**12) - 2*(sigma**6/r**6))

    return E_v

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
