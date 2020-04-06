
import os
import hoomd
import hoomd.md
import math
import random
import numpy
import gsd
import gsd.hoomd
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
circles = False


def hardContact(r, rmin, ramx, l_0,K):
    V = K*(r-l_0)**2
    F = -2*K*(r-l_0)
    return (V,F)

def harmonicF(r, rmin, rmax,l_0,l_max,K):
    if r < (l_0-l_max + 0.00001):
        V = K*(r-l_0)**2
        F = -2*K*(r-l_0)
        return (V,F)
    if r >= (l_0 + l_max):
        V = - 0.5*K*l_max**2*math.log(1-(l_max-0.001)**2/l_max**2)
        F =  - (K*(l_max-0.001))/(1-(l_max-0.001)**2/l_max**2)
        # V = K*(r-l_0)**6
        # F = -5*K*(r-l_0)**5
        return (V,F)

    V = - 0.5*K*l_max**2*math.log(1-(r-l_0)**2/l_max**2)
    F =  - (K*(r-l_0))/(1-(r-l_0)**2/l_max**2)
    return (V,F)


class PolymerSimulationParameters():
    def __init__(self,data = None):
        if data == None:
            print("********\n*\n*NO PARAMETERS HAVE BEEN PASSED INTO POLYMER_SIMULATION_PARAMETERS\n*initializing all parameters to 0\n*\n********")
            self.sheerForce=0
            self.length=int(0)
            self.lines=int(0)
            self.rez=0
            self.K=0
            self.l_0=0
            self.l_max=0
            self.pull=0
            self.amplitude=0
            self.gamma=0
            self.kbT=0
            self.dt=0
            self.probePeriod=int(0)
            self.runLength=int(0)

        else:
            a=0
            self.sheerForce = data[a]
            a+=1
            self.length     = int(data[a])
            a+=1
            self.lines      = int(data[a])
            a+=1
            self.rez        = data[a]
            a+=1
            self.K          = data[a]
            a+=1
            self.l_0        = data[a]
            a+=1
            self.l_max      = data[a]
            a+=1
            self.pull       = data[a]
            a+=1
            self.amplitude  = data[a]
            a+=1
            self.gamma      = data[a]
            a+=1
            self.kbT        = data[a]
            a+=1
            self.dt         = data[a]
            a+=1
            self.probePeriod= int(data[a])
            a+=1
            self.runLength  = int(data[a])

    def setSheerForce(self,x):
        self.sheerForce = x
    def setLength(self,x):
        self.length = int(x)
    def setNumberChains(self,x):
        self.lines = int(x)
    def setPairRadius(self,x):
        self.rez = x
    def setPairPotentialStrength(self,x):
        self.K = x
    def setPairRadiusEqualibrium(self,x):
        self.l_0 = x
    def setPairMaximumRadius(self,x):
        self.l_max = x
    def setForcePull(self,x):
        self.pull = x
    def setPeriodicAmplitude(self,x):
        self.amplitude = x
    def setGamma(self,x):
        self.gamma = x
    def setKBT(self,x):
        self.kbT = x
    def setTimeStep(self,x):
        self.dt = x
    def setProbePeriod(self,x):
        self.probePeriod = int(x)
    def setRunLength(self,x):
        self.runLength = int(x)


    def getSheerForce(self):
        return self.sheerForce
    def getLength(self):
        return self.length
    def getNumberChains(self):
        return self.lines
    def getPairRadius(self):
        return self.rez
    def getPairPotentialStrength(self):
        return self.K
    def getPairRadiusEqualibrium(self):
        return self.l_0
    def getPairMaximumRadius(self):
        return self.l_max
    def getPullForce(self):
        return self.pull
    def getPeriodicAmplitude(self):
        return self.amplitude
    def getGamma(self):
        return self.gamma
    def getKBT(self):
        return self.kbT
    def  getTimeStep(self):
        return self.dt
    def  getProbePeriod(self):
        return self.probePeriod
    def  getRunLength(self):
        return self.runLength

class PolymerSimulation():
    def __init__(self):
        self.parameter = None
        random.seed()

    def init(self, name=None, parameter=None,initializer='--mode=cpu'):
        hoomd.context.initialize(initializer)
        self.parameter = parameter
        random.seed()
        if name is None:

            self.setupFileSystem(True)
            #1. create new directory for this run

            self.populateSystem()
            #2. load and populate snapshot for simulation

            self.initializeForces()
            self.initializeIntegrator()

            #3. insert parameters for simulation
        else:
            pass
            ## TODO: Add ability to load snapshot

    def run(self,forceRange=None,equalibriate = False):
        self.setupFileSystem()
        for i in range(forceRange[1]):
            

    def initializeIntegrator(self):

        hoomd.md.integrate.mode_standard(dt=self.parameter.getTimeStep());
        bs = hoomd.md.integrate.brownian(group=self.all, kT=self.parameter.getKBT(), seed=random.randint(0,999999));
        bs.set_gamma('A', gamma=self.parameter.getGamma())
        bs.set_gamma('B', gamma=self.parameter.getGamma())
        bs.set_gamma('C', gamma=self.parameter.getGamma())


    def setupFileSystem(self,init=False,name=None):
        from datetime import datetime
        time = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')

        if init is True:
            if name==None:
                self.DirectoryName = "Simulation_" + time
                os.system("mkdir " + self.DirectoryName)
            else:
                self.DirectoryName = "Simulation_" + name
                os.system("mkdir " + self.DirectoryName)


        else:
            tempName = "Simulation_" + time
            self.currentSimulation = tempName
            os.system("mkdir " + self.DirectoryName + "/" + tempName)



    def initializeForces(self):
        added = 0
        lines       = self.parameter.getNumberChains()
        length      = self.parameter.getLength()
        l_0         = self.parameter.getPairRadiusEqualibrium()
        l_max       = self.parameter.getPairMaximumRadius()
        K           = self.parameter.getPairPotentialStrength()
        amplitude   = self.parameter.getPeriodicAmplitude()
        pull        = self.parameter.getPullForce()
        rmax        = l_0 + l_max - 0.000000000001

        self.all    = hoomd.group.all()
        self.anchor = hoomd.group.type(name='anchor', type='A')
        self.pulley = hoomd.group.type(name='pulley',type='B')
        self.chain  = hoomd.group.type(name='chain', type='C')

        nl = hoomd.md.nlist.cell();



        table = hoomd.md.pair.table(width=10000, nlist=nl) #I think the radius cut is the ran
        table.pair_coeff.set('A', 'A', func=hardContact,rmin=0,rmax=l_0, coeff=dict(l_0=l_0, K=K))
        table.pair_coeff.set('A', 'B', func=hardContact,rmin=0,rmax=l_0, coeff=dict(l_0=l_0, K=K))
        table.pair_coeff.set('B', 'B', func=hardContact,rmin=0,rmax=l_0, coeff=dict(l_0=l_0, K=K))

        table.pair_coeff.set('B', 'C', func=hardContact,rmin=0,rmax=l_0, coeff=dict(l_0=l_0, K=K))
        table.pair_coeff.set('A', 'C', func=hardContact,rmin=0,rmax=l_0, coeff=dict(l_0=l_0, K=K))
        table.pair_coeff.set('C', 'C', func=hardContact,rmin=0,rmax=l_0, coeff=dict(l_0=l_0, K=K))

        nl.reset_exclusions(exclusions = [])

        harmonic = hoomd.md.bond.table(width = 100000);
        harmonic.bond_coeff.set('polymer', func = harmonicF,rmin=0, rmax=10,coeff = dict(l_0=l_0,l_max=l_max,K=K));

        # fene = hoomd.md.bond.fene()
        # fene.bond_coeff.set('polymer', k=10**3, r0=6.0, sigma=0.0, epsilon= 0.0)

        self.tensionForce = hoomd.md.force.constant(group = self.pulley, fvec=(0.0,pull,0.0))  # ?
        self.sheer = hoomd.md.force.constant(group = self.anchor, fvec=(0.0,0.0,0.0))
        periodic = hoomd.md.external.periodic()
        periodic.force_coeff.set('A', A=amplitude, i=0, w=1, p=lines+added)
        periodic.force_coeff.set('B', A=amplitude, i=0, w=1, p=lines+added)
        periodic.force_coeff.set('C', A=amplitude, i=0, w=1, p=lines+added)

        periodic.force_coeff.set('A', A=-10000000.0, i=1, w=1, p=100)

    def definePositions(self):
        lines   = self.parameter.getNumberChains()
        length  = self.parameter.getLength()
        l_0     = self.parameter.getPairRadiusEqualibrium()
        l_max   = self.parameter.getPairMaximumRadius()
        pos     = []

        for i in range(lines):
            x = i - lines/2
            y = 0

            for j in range(length):
                y = y + l_0
                pos.append([x,y,0])

        return pos

    def defineBonds(self):
        lines   = self.parameter.getNumberChains()
        length  = self.parameter.getLength()
        bonds   = []

        for i in range(lines):
            for j in range(1,length):
                bonds.append([(j-1 + i*length),(j + i*length)])
        return bonds


    def defineIds(self):
        lines   = self.parameter.getNumberChains()
        length  = self.parameter.getLength()
        id     = []

        for i in range(lines):
            for j in range(length):
                if(j==0):
                    id.append(0)
                elif (j==(length-1)):
                    id.append(1)
                else:
                    id.append(2)
        return id

    def defineMass(self):
        mass    = []
        lines   = self.parameter.getNumberChains()
        length  = self.parameter.getLength()

        for i in range(lines*length):
            mass.append(1)
        return mass


    def populateSystem(self):
        pos = []
        bonds = []
        id = []
        mass = []
        types = ['A','B','C']

        self.boxdim=[0,0]
        self.boxdim[0] = self.parameter.getNumberChains()
        self.boxdim[1] = self.parameter.getLength()*16

        boundary = hoomd.data.boxdim(Lx = self.boxdim[0], Ly = self.boxdim[1], dimensions=2)
        snapshot = hoomd.data.make_snapshot(N=  int(self.parameter.getNumberChains()  *  self.parameter.getLength()), box=boundary, particle_types=types, bond_types=['polymer'])

        pos     = self.definePositions()
        bonds   = self.defineBonds()
        id      = self.defineIds()
        mass    = self.defineMass()


        snapshot.particles.position[:] = pos
        snapshot.bonds.resize(len(bonds))
        snapshot.bonds.group[:] = bonds
        snapshot.particles.typeid[:] = id
        snapshot.particles.mass[:] = mass

        hoomd.init.read_snapshot(snapshot)

    def simulationReadMeDump(self,force = None,name=None, dir=None):
        text = None
        if name != None and dir !=None:
            text = open(dir + name + "_simulation_parameters.txt","w+")
            text.write("sheerForce=" + str(force) + "\n")
        else:
            text = open(self.DirectoryName + "/" + self.currentSimulation + "/simulation_parameters.txt","w+");
        text.write("length="    +      str(self.parameter.getLength())                 + "\n")
        text.write("lines="     +      str(self.parameter.getNumberChains())           + "\n")
        text.write("rez="       +      str(self.parameter.getPairRadius())             + "\n")
        text.write("K="         +      str(self.parameter.getPairPotentialStrength())  + "\n")
        text.write("l_0="       +      str(self.parameter.getPairRadiusEqualibrium())  + "\n")
        text.write("l_max="     +      str(self.parameter.getPairMaximumRadius())      + "\n")
        text.write("pull="      +      str(self.parameter.getPullForce())              + "\n")
        text.write("amplitude=" +      str(self.parameter.getPeriodicAmplitude())      + "\n")
        text.write("gamma="     +      str(self.parameter.getGamma())                  + "\n")
        text.write("kbT="       +      str(self.parameter.getKBT())                    + "\n")
        text.write("dt="        +      str(self.parameter.getTimeStep())               + "\n")
        text.write("probePeriod="+     str(self.parameter.getProbePeriod())            + "\n")
        text.write("runLength=" +      str(self.parameter.getRunLength())              + "\n")
        text.write("BoxDim="    +      str(self.boxdim[0]) + "," + str(self.boxdim[1])  + "\n")
        text.close()



class DataVisualizer():
        def __init__(self,basedirectory="",interval=0):

            self.interval = interval
            self.directory = basedirectory


        def init(self):
            import glob

            location = self.directory
            gsdFileLocations = [file for file in glob.glob(location + "**/*.gsd", recursive=True)]
            simulationParameterLocations = [file for file in glob.glob(location + "**/*.txt", recursive=True)]
            energyFileLocations = [file for file in glob.glob(location + "**/*.log", recursive=True)]

            for i in range(len(simulationParameterLocations)):

                self.name = str(gsdFileLocations[i]).split('/')[-1].split('.')[0]
                location = str(gsdFileLocations).split('/')[0] + "/" + str(gsdFileLocations).split('/')[1] + "/" + str(gsdFileLocations).split('/')[2] + "/ "
                self.getSimulationParameters(simulationParameterLocations[i])
                self.constructPolymerObjects(gsdFileLocations[i])
                #self.plotGeneralPolymerProfiles(saveLocation = location)
                #self.getPositionProbabilityData(name = self.name,location=location,Interval = self.interval)
                #self.plotEnergy(location=energyFileLocations[i],name=self.name+"energy.png")
                #self.animatePositionData(fps=300)

        def plotEnergy(self,location="",name='energy.png'):
            print(location)
            data = numpy.genfromtxt(fname=location, skip_header=True);
            plt.figure();
            plt.clf()
            plt.plot(data[:,0], data[:,1]);
            plt.xlabel('time step');
            plt.ylabel('potential_energy');
            plt.savefig(name)

        def generatePotential(self,dimx,rez):
            boxdimx = dimx
            rez=rez
            p=[]
            for i in range(rez):
                p.append([])
                for j in range(rez):
                    p[-1].append(-math.cos(j*2*math.pi/rez*boxdimx - boxdimx/2))

            return p


        def animatePositionData(self,fps=500,Interval = 0):
            circles = False
            fig, ax = plt.subplots(figsize=(10, 10))
            xmax = self.gsd_data[0].particles.position[:,0][0]
            xmin = xmax
            ymax = self.gsd_data[0].particles.position[:,1][0]
            ymin = ymax

            for i in range(len(self.gsd_data)):
                for j in range(len(self.gsd_data[int(i)].particles.position[:,0])):
                    if self.gsd_data[int(i)].particles.position[:,0][j] < xmin:
                        xmin = self.gsd_data[int(i)].particles.position[:,0][j]
                    if self.gsd_data[int(i)].particles.position[:,0][j] > xmax:
                        xmax = self.gsd_data[int(i)].particles.position[:,0][j]

                    if self.gsd_data[int(i)].particles.position[:,1][j] < ymin:
                        ymin = self.gsd_data[int(i)].particles.position[:,1][j]
                    if self.gsd_data[int(i)].particles.position[:,1][j] > ymax:
                        ymax = self.gsd_data[int(i)].particles.position[:,1][j]
            potential = self.generatePotential(self.parameters.getNumberChains(),1000)
            print("Total frames = " + str(len(self.gsd_data)))
            def update(i):

                # Update the line and the axes (with a new xlabel). Return a tuple of
                # "artists" that have to be redrawn for this frame.
                print(int(i))
                x = self.gsd_data[int(i)].particles.position[:,0]
                y = self.gsd_data[int(i)].particles.position[:,1]

                ax.clear()
                #ax.set_size_inches(10,10)
                ax.imshow(potential,extent=[-self.parameters.getNumberChains()/2,self.parameters.getNumberChains()/2,ymin,ymax])
                for i in range(self.parameters.getNumberChains()):
                    length = self.parameters.getLength()
                    ax.plot(x[i*length:i*length + length],y[i*length:i*length + length])
                    ax.plot(x[i*length:i*length + length],y[i*length:i*length + length],'o')

                ax.set_xlim(xmin,xmax)
                ax.set_ylim(ymin,ymax)
                ax.set_aspect('auto')
                ax.set_xlabel(str(i) + "/" + str(len(self.gsd_data)))

                return ax
            anim = FuncAnimation(fig, update, frames=numpy.arange(int(len(self.gsd_data)*Interval), len(self.gsd_data)), interval=int(fps))
            anim.save(self.name +'.gif', dpi=80, writer='imagemagick')
            print("finished")


        def getSimulationParameters(self,fileLocation):
            file = open(fileLocation,'r');
            data = []
            lines = file.readlines()
            for i in range(len(lines)-1):
                data.append(float(lines[i].split('=')[1]))
            data.append([float(lines[len(lines)-1].split('=')[1].split(',')[0]),float(lines[len(lines)-1].split('=')[1].split(',')[1])])

            self.parameters = PolymerSimulationParameters(data=data)

            file.close()

        def constructPolymerObjects(self,fileLocation):

            self.gsd_data =gsd.hoomd.open(fileLocation,'rb')
            self.polymers = []
            print(len(self.gsd_data))
            #print(self.parameters.getNumberChains(),self.parameters.getLength(),len(self.gsd_data))
            for i in range(int(self.parameters.getNumberChains())):

                x = []
                y = []
                for j in range(int(self.parameters.getLength())):

                    x.append([])
                    y.append([])
                    for k in range(int(len(self.gsd_data)*self.interval),len(self.gsd_data)):

                        x[-1].append(self.gsd_data[k].particles.position[i*self.parameters.getLength() + j,0])

                        y[-1].append(self.gsd_data[k].particles.position[i*self.parameters.getLength() + j,1])
                    #print(i*self.parameters.getNumberChains()+j)
                self.polymers.append(PolymerObject([x,y],L=self.parameters.getNumberChains()))


        def gaussian(self,x, mu, sig):
            return numpy.exp(-numpy.power(x - mu, 2.) / (2 * numpy.power(sig, 2.)))

        def plotGeneralPolymerProfiles(self,savePlot = True,saveLocation=""):

            if savePlot == True:
                for i in range(len(self.polymers)):
                    plt.clf()
                    print(i)
                    generalProfileWidths = self.polymers[i].getWidths()
                    #y = range(len(generalProfileWidths),0,-1)
                    y=[]
                    x=[]
                    mu = self.polymers[i].getMean()
                    for j in range(len(generalProfileWidths)):
                        y.append(mu[j][1])
                        x.append(mu[j][0])

                    xerr = []
                    for j in range(len(y)):
                        xerr.append(generalProfileWidths[j][0])
                    plt.plot(x,y)
                    plt.errorbar(x,y,xerr=xerr,fmt='o')
                    #t = numpy.linspace(start=x[0]-1, stop=x[-1]+1, num=50)
                    for j in range(len(y)):
                        pass
                        #plt.plot(t, self.gaussian(t, x[j] , generalProfileWidths[j][0]) + y[j],color='k',linewidth = 1)
                    print(self.name + "_" + str(i)+"_widths.png")
                    plt.savefig(self.name + "_" + str(i)+"_widths.png")
                    #os.system("mv " + self.name + "_p" + str(i)+ ".png " + saveLocation)






        def getPositionProbabilityData(self,rez=None,Interval=0,name = "foo",location=""):
            if rez == None:
                rez = [100,100]


            p = []

            p_t = []
            p_t.append([])
            p_t.append([])

            for particle in range(len(self.gsd_data[0].particles.position)):
                p.append([])
                p[-1].append([])
                p[-1].append([])
                for i in range(int(len(self.gsd_data)*Interval),len(self.gsd_data)):
                    p_t[0].append(self.gsd_data[i].particles.position[particle,0])
                    p_t[1].append(self.gsd_data[i].particles.position[particle,1])



            #Plot heatmap
            plt.clf()

            plt.title('heatmap test')
            plt.ylabel('y')
            plt.xlabel('x')
            plt.hist2d(p_t[0],p_t[1],bins=(rez[0],rez[1]))
            writename = "heatmap_" +name
            plt.savefig(writename + ".png")
            os.system("mv " +name+".png " + location)

class GlobalDataAnalyzer():
    def __init__(self,location):

        import glob

        gsdFileLocations = [file for file in glob.glob(location + "**/*.gsd", recursive=True)]
        simulationParameterLocations = [file for file in glob.glob(location + "**/*.txt", recursive=True)]


        self.parameters = []
        self.polymers = []
        self.interval = 0.5
        for i in range(len(simulationParameterLocations)):

            self.name = str(gsdFileLocations[i]).split('/')[-1].split('.')[0]
            location = str(gsdFileLocations).split('/')[0] + "/" + str(gsdFileLocations).split('/')[1] + "/" + str(gsdFileLocations).split('/')[2] + "/" + str(gsdFileLocations).split('/')[3] + "/"
            print(location)
            self.getSimulationParameters(simulationParameterLocations[i])
            self.constructPolymerObjects(gsdFileLocations[i])


    def getSimulationParameters(self,fileLocation):
        file = open(fileLocation,'r');
        data = []
        lines = file.readlines()
        print(lines)
        for i in range(len(lines)-1):
            data.append(float(lines[i].split('=')[1]))
        data.append([float(lines[len(lines)-1].split('=')[1].split(',')[0]),float(lines[len(lines)-1].split('=')[1].split(',')[1])])

        self.parameters.append(PolymerSimulationParameters(data=data))

        file.close()

    def plotTiltbyForce(self,tension, temp):
        tilt = []
        x = []
        sheerForce = []
        for i in range(len(self.polymers)):
            if (self.parameters[i].getPullForce() == tension and self.parameters[i].getKBT() == temp):
                sheerForce.append(self.parameters[i].getSheerForce())

                tilt.append(0)
                x.append(0)
                for j in range(len(self.polymers[i])):
                    tilt[-1] += self.polymers[i][j].getTilt()
                    x[-1] += self.polymers[i][j].getDX()
                tilt[-1] = tilt[-1]/len(self.polymers[i])
                #print(x[-1])
                x[-1] = x[-1]/len(self.polymers[i])

        plt.clf()

        plt.title('Tilt vs SheerForce: Ft:' + str(tension) + ' kbT=' + str(temp))
        #plt.plot(sheerForce,tilt,'.')
        plt.plot(sheerForce,x,'.')
        plt.savefig("forceVstilt.png")




    def constructPolymerObjects(self,fileLocation):

        self.gsd_data =gsd.hoomd.open(fileLocation,'rb')
        self.polymers.append([])

        #print(self.parameters.getNumberChains(),self.parameters.getLength(),len(self.gsd_data))
        for i in range(int(self.parameters[-1].getNumberChains())):

            x = []
            y = []
            for j in range(int(self.parameters[-1].getLength())):

                x.append([])
                y.append([])
                for k in range(int(len(self.gsd_data)*self.interval),len(self.gsd_data)):

                    x[-1].append(self.gsd_data[k].particles.position[i*self.parameters[-1].getLength() + j,0])

                    y[-1].append(self.gsd_data[k].particles.position[i*self.parameters[-1].getLength() + j,1])
                #print(i*self.parameters[-1].getNumberChains()+j)
            self.polymers[-1].append(PolymerObject([x,y],L=self.parameters[-1].getNumberChains()))




class PolymerObject():
    def __init__(self,data,L = 0):
        self.L = L
        self.particle = data
        self.dx=0
        # x = []
        # y = []
        #
        # for i in range(len(data[0])):
        #         x.append(data[0][i][10])
        #         y.append(data[1][i][10])
        # plt.plot(x,y,'o')
        # plt.plot(x,y)
        # plt.savefig("TEST.png")
        self.generateWidths()
    def getParticleData(self,id):
        return [self.particle[0][id],self.particle[1][id]]

    def getParticleWidth(self,id):
        return self.particleWidth[id]

    def generateWidths(self):
        self.particleWidth = []
        self.particleMean = []
        self.base = self.particle[0][0][0]
        for i in range(len(self.particle[0])):
            self.particleWidth.append([   self.stdDevX(self.particle[0][i])  , self.stdDevY(self.particle[1][i])])
            self.particleMean.append([self.calculateMeanX(self.particle[0][i]),self.calculateMeanY(self.particle[1][i])])
        for i in range(1,len(self.particleMean)):
            if self.particleMean[i][0] - self.particleMean[i-1][0] <= -self.L/2:
                self.particleMean[i][0] = self.particleMean[i][0] + self.L
            elif self.particleMean[i][0] - self.particleMean[i-1][0] >= self.L/2:
                self.particleMean[i][0] = self.particleMean[i][0] - self.L

        self.tilt = numpy.arctan((self.particleMean[-1][1] - self.particleMean[0][1])/(self.particleMean[-1][0] - self.particleMean[0][0]))/math.pi*180
        self.dx = (self.particleMean[-1][0] - self.particleMean[0][0])
        print(self.dx)
        self.Length = math.sqrt((self.particleMean[-1][1] - self.particleMean[0][1])**2 + (self.particleMean[-1][0] - self.particleMean[0][0])**2)

    def getLength(self):
        return self.Length
    def calculateMeanX(self,arr):
        for i in range(1,len(arr)):
            if arr[i] - arr[i-1] <= -self.L/2:
                arr[i] = arr[i] + self.L
            elif arr[i] - arr[i-1] >= self.L/2:
                arr[i] = arr[i] - self.L

        return sum(arr)/len(arr)
    def getDX(self):
        return self.dx
    def getTilt(self):
        return self.tilt
    def calculateMeanY(self,arr):
        return sum(arr)/len(arr)

    def getMean(self):
        return self.particleMean
    def getWidths(self):
        return self.particleWidth

    def stdDevX(self,data):
        std = 0
        mean = self.calculateMeanX(data)
        difference = 0
        for i in range(1,len(data)):
            if data[i] - data[i-1] < -self.L/2:
                data[i] = data[i] + self.L
            elif data[i] - data[i-1] > self.L/2:
                data[i] = data[i] - self.L
        for i in range(len(data)):
            difference += (data[i]-mean)**2
        std = math.sqrt(1/(len(data)-1)*difference)
        return std
    def stdDevY(self,data):
        std = 0
        mean = sum(data)/len(data)
        difference = 0
        for i in range(len(data)):

            difference += (data[i]-mean)**2
        std = math.sqrt(1/(len(data)-1)*difference)

        return std
