import hoomd
from hoomd import md
import hoomd.md.external
import math
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
import gsd
import gsd.hoomd
import sys
import numpy
import seaborn

class Analyze():
    """docstring for Analyse."""

    def __init__(self, gsd_directory="",names=None):
        if(names == None):
            self.gsd_directory = gsd_directory
            self.gsd_data = gsd.hoomd.open(self.gsd_directory + "polymer.gsd", 'rb');
            self.energy_data = numpy.genfromtxt(fname="energy.log", skip_header=True);
            self.simulation_data = self.getSimulationParameters("simulation_parameters.txt")
        else:
            self.gsd_directory = gsd_directory
            self.gsd_data = gsd.hoomd.open(self.gsd_directory  + names[0] + ".gsd", 'rb');
            #self.energy_data = numpy.genfromtxt(fname=names[1]+".log", skip_header=True);
            self.simulation_data = self.getSimulationParameters("simulation_parameters.txt")


    def getPositionProbabilityData(self,rez=None,indvidual=None,Interval=0,name = "foo.pdf"):
        if rez == None:
            rez = [int(self.boxdim[0]),self.boxdim[1]]
        from matplotlib.backends.backend_pdf import PdfPages
        #e_ave = sum(self.energy_data[:,1])/len(self.energy_data[:,1])
        #dt = 2*e_ave
        px = []
        py = []
        for i in range(len(self.gsd_data)):
            px.append(self.gsd_data[i].particles.position[:,0])
            py.append(self.gsd_data[i].particles.position[:,1])

        p = []

        p_t = []
        p_t.append([])
        p_t.append([])

        for particle in range(len(px[0])):
            #ERROR HERE
            p.append([])
            p[-1].append([])
            p[-1].append([])
            for i in range(int(len(self.gsd_data)*Interval),len(self.gsd_data)):
                p_t[0].append(self.gsd_data[i].particles.position[particle,0])
                p_t[1].append(self.gsd_data[i].particles.position[particle,1])
                if indvidual!= None:
                    p[-1][0].append(self.gsd_data[i].particles.position[particle,0])
                    p[-1][1].append(self.gsd_data[i].particles.position[particle,1])
        with PdfPages(name + ".pdf") as pdf:
            if indvidual!= None:
                for particle in range(len(px[0])):
                    ax = plt.figure(particle)
                    heatmap, xedges, yedges = numpy.histogram2d(p[particle][0], p[particle][1], bins=(rez[0],rez[1]))
                    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
                    # Plot heatmap
                    plt.clf()
                    plt.title(' heatmap test')
                    plt.ylabel('y')
                    plt.xlabel('x')
                    plt.imshow(heatmap,aspect='auto',extent=extent)
                    pdf.savefig(ax)
                    plt.close(ax)
            ax = plt.figure(len(px[0]))
            heatmap, xedges, yedges = numpy.histogram2d(p_t[1], p_t[0], bins=(rez[0],rez[1]))
            extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

            #Plot heatmap
            plt.clf()
            plt.title('heatmap test')
            plt.ylabel('y')
            plt.xlabel('x')
            plt.hist2d(p_t[0],p_t[1],bins=(200,200))
            pdf.savefig(ax)





    def getSimulationParameters(self, filename):
        file = open(filename);

        lines = file.readlines()

        self.name = lines[0].split(':')[1]
        self.period = int(lines[1].split(":")[1])
        self.einterval = int(lines[2].split(":")[1])
        self.dt = float(lines[3].split(":")[1])
        self.runLength = int(lines[4].split(":")[1])
        self.boxdim = [int(lines[5].split(":")[1].split(',')[0]),int(lines[5].split(":")[1].split(',')[1])]

    def plotPositionData(self,fps=500,Interval = 0):
        fig, ax = plt.subplots()
        fig.set_tight_layout(True)
        a,b = self.find_max(self.gsd_data)
        ax.set_xlim([-a,a])
        ax.set_ylim([0,b])

        def update(i):

            # Update the line and the axes (with a new xlabel). Return a tuple of
            # "artists" that have to be redrawn for this frame.

            x = self.gsd_data[int(i)].particles.position[:,0]
            y = self.gsd_data[int(i)].particles.position[:,1]
            ax.clear()
            ax.set_xlim(-a,a)
            ax.set_ylim(-b,b)
            ax.plot(x,y)
            ax.plot(x,y,'o')
            ax.set_xlabel(str(i) + "/" + str(len(self.gsd_data)))
            return ax

        anim = FuncAnimation(fig, update, frames=numpy.arange(int(len(self.gsd_data)*Interval), len(self.gsd_data)), interval=int(fps))
        anim.save('Polymer.gif', dpi=80, writer='imagemagick')
        print("finished")






    def find_max(self,arr):
        xm = 0
        ym = 0
        for i in range(len(arr)):
            if max(abs(arr[i].particles.position[:,0]) > xm):
                xm = max(abs(arr[i].particles.position[:,0]))
            if max(abs(arr[i].particles.position[:,1]) > ym):
                ym = max(abs(arr[i].particles.position[:,1]))
        return xm, ym
