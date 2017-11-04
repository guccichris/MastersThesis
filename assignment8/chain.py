# File Name: chain.py
# Description: A chain of bonded atoms floating above a rectangular substrate.
# Author: Christopher Parker
# Created: Fri Nov 03, 2017 | 10:32P EDT
# Last Modified: Fri Nov 03, 2017 | 11:00P EDT

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
from two_atom_spring import spring_and_VDWForce


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
fleaters = np.zeros(3)

for i in range(len(floaters)):
    floaters[i] = (i, i, 1)


