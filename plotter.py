
import hoomd
from hoomd import md
import hoomd.md.external
import math
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
import gsd
import gsd.hoomd
import sys
import numpy as np

def find_max(arr):
    xm = 0
    ym = 0
    for i in range(len(arr)):
        if max(abs(arr[i].particles.position[:,0]) > xm):
            xm = max(abs(arr[i].particles.position[:,0]))
        if max(abs(arr[i].particles.position[:,1]) > ym):
            ym = max(abs(arr[i].particles.position[:,1]))
    return xm, ym



t = gsd.hoomd.open(sys.argv[1] + ".gsd", 'rb');

fig, ax = plt.subplots()
fig.set_tight_layout(True)
a,b = find_max(t)
ax.set_xlim([-a,a])
ax.set_ylim([0,b])

outF = open("distance-debug.txt","w")
for a in range(len(t)):
    outF.write("\n\n\n\n\n")
    for j in range(len(t[a].particles.position)-1):
        x = t[a].particles.position[j,1]
        x1 = t[a].particles.position[j+1,1]
        outF.write(str(x1-x))
outF.close()



def update(i):

    # Update the line and the axes (with a new xlabel). Return a tuple of
    # "artists" that have to be redrawn for this frame.

    x = t[int(i)].particles.position[:,0]
    y = t[int(i)].particles.position[:,1]
    ax.clear()
    ax.set_xlim(-a,a)
    ax.set_ylim(-b,b)
    ax.plot(x,y)
    ax.plot(x,y,'o')
    ax.set_xlabel(str(i) + "/" + str(len(t)))
    return ax



if (sys.argv[2]=="-a"):
    anim = FuncAnimation(fig, update, frames=np.arange(0, len(t)), interval=int(sys.argv[3]))
    anim.save('Polymer.gif', dpi=80, writer='imagemagick')

    print("finished")





else:

    print(len(t))
    print(len(t[3].particles.position))
    print(t[3].particles.position)



    for i in range(len(t)):
        x = t[i].particles.position[:,0]
        y = t[i].particles.position[:,1]
        plt.figure(i)
        plt.plot(x,y)
        plt.show()


# import numpy
# data = numpy.genfromtxt(fname='log-output.log', skip_header=True);
#
#
# plt.figure(figsize=(4,2.2), dpi=140);
# plt.plot(data[:,0], data[:,1]);
# plt.xlabel('time step');
# plt.ylabel('potential_energy [kJ/mol]');
#
#
#
# plt.figure(figsize=(4,2.2), dpi=140);
# plt.plot(data[:,0], data[:,2]);
# plt.xlabel('time step');
# plt.ylabel('temperature');
#
# plt.show()
