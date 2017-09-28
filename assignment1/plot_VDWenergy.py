import matplotlib.pyplot as plt
import numpy as np
from truncatedSum_Ev import VDW


floater = (np.linspace(0,10,500),1,.8)
V = []
if type(floater[0]).__module__ == np.__name__:
    for pt in floater[0]:
        V.append(VDW(pt,floater[1],floater[2]))
    plt.plot(floater[0],V,'k')
    plt.xlabel('x (y = 1, z = 1)')
    plt.ylabel('VDW energy')
    plt.savefig('change_in_x.png')
    plt.show()
elif type(floater[1]).__module__ == np.__name__:
    for pt in floater[1]:
        V.append(VDW(floater[0],pt,floater[2]))
    plt.plot(floater[1],V,'k')
    plt.xlabel('y (x = 1, z = 1)')
    plt.ylabel('VDW energy')
    plt.savefig('change_in_y.png')
    plt.show()
elif type(floater[2]).__module__ == np.__name__:
    for pt in floater[2]:
        V.append(VDW(floater[0],floater[1],pt))
    plt.plot(floater[2],V,'k')
    plt.xlabel('z (x = 1, y = 1)')
    plt.ylabel('VDW energy')
    plt.savefig('change_in_z.png')
    plt.show()
else:
    print("Only one point given!")

# try setting only one of x, y or z and letting the other two be linspcaes. then plot in 3d
# mayavi (try looking into this for 3d plotting)
