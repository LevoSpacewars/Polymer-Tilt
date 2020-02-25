
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
    def __init__(self,name=None,parameter=None,initializer='--mode=gpu'):
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

    def run(self,forceRange=None):
        self.setupFileSystem()
        for i in range(int(forceRange[2])):
            conv = round((forceRange[1]-forceRange[0])/forceRange[2]*i + forceRange[0],int(math.log(forceRange[2])))
            gsdname="polymer_" + str(conv).replace(".","_") + ".gsd"
            energyname = "energy_" + str(conv).replace(".","_") + ".log"
            hoomd.analyze.log(filename=energyname,
                          quantities=['potential_energy', 'temperature'],
                          period=self.parameter.getProbePeriod(),
                          overwrite=True);

            self.tensionForce = hoomd.md.force.constant(group=self.pulley , fvec=(conv,self.parameter.getPullForce(),0.0))

            self.parameter.setSheerForce(conv)

            hoomd.dump.gsd(gsdname, period=self.parameter.getProbePeriod(), group=self.all, overwrite=True);
            hoomd.run(self.parameter.getRunLength())


            dirname = "simulationForce_" + str(conv).replace(".","_")
            os.system("mkdir " + self.DirectoryName + "/" + self.currentSimulation + "/" +dirname)
            #print(self.DirectoryName + "/" + self.currentSimulation + "/" + dirname + "/")
            self.simulationReadMeDump(force = conv, name = gsdname[0:-4], dir = self.DirectoryName + "/" + self.currentSimulation + "/" + dirname + "/")

            os.system("mv " + gsdname + " " + self.DirectoryName + "/" + self.currentSimulation + "/" + dirname + "/")
            os.system("mv " + energyname+ " " +self.DirectoryName + "/" + self.currentSimulation + "/" + dirname + "/")
        return self.DirectoryName + "/"

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
        added = 1
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



        table = hoomd.md.pair.table(width=1000000, nlist=nl) #I think the radius cut is the ran
        table.pair_coeff.set('A', 'A', func=hardContact,rmin=0,rmax=l_0, coeff=dict(l_0=l_0, K=K))
        table.pair_coeff.set('A', 'B', func=hardContact,rmin=0,rmax=l_0, coeff=dict(l_0=l_0, K=K))
        table.pair_coeff.set('B', 'B', func=hardContact,rmin=0,rmax=l_0, coeff=dict(l_0=l_0, K=K))

        table.pair_coeff.set('B', 'C', func=hardContact,rmin=0,rmax=l_0, coeff=dict(l_0=l_0, K=K))
        table.pair_coeff.set('A', 'C', func=hardContact,rmin=0,rmax=l_0, coeff=dict(l_0=l_0, K=K))
        table.pair_coeff.set('C', 'C', func=hardContact,rmin=0,rmax=l_0, coeff=dict(l_0=l_0, K=K))

        nl.reset_exclusions(exclusions = [])

        harmonic = hoomd.md.bond.table(width = 10000000);
        harmonic.bond_coeff.set('polymer', func = harmonicF,rmin=0, rmax=10,coeff = dict(l_0=l_0,l_max=l_max,K=K));

        # fene = hoomd.md.bond.fene()
        # fene.bond_coeff.set('polymer', k=10**3, r0=6.0, sigma=0.0, epsilon= 0.0)

        self.tensionForce = hoomd.md.force.constant(group = self.pulley, fvec=(0.0,pull,0.0))  # ?
        periodic = hoomd.md.external.periodic()
        periodic.force_coeff.set('A', A=amplitude, i=0, w=1, p=lines+added)
        periodic.force_coeff.set('B', A=amplitude, i=0, w=1, p=lines+added)
        periodic.force_coeff.set('C', A=amplitude, i=0, w=1, p=lines+added)

        periodic.force_coeff.set('A', A=-10000000.0, i=1, w=1, p=5)


    def definePositions(self):
        lines   = self.parameter.getNumberChains()
        length  = self.parameter.getLength()
        l_0     = self.parameter.getPairRadiusEqualibrium()
        l_max   = self.parameter.getPairMaximumRadius()
        pos     = []

        for i in range(lines):
            x = random.uniform(0, lines) - lines/2
            y = 0

            for j in range(length):
                y = y + random.uniform(l_0, l_max)
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

        boundary = hoomd.data.boxdim(Lx = self.boxdim[0]+1, Ly = self.boxdim[1], dimensions=2)
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
            import glob

            location = basedirectory
            gsdFileLocations = [file for file in glob.glob(location + "**/*.gsd", recursive=True)]
            simulationParameterLocations = [file for file in glob.glob(location + "**/*.txt", recursive=True)]
            #print(simulationParameterLocations)
            #print(len(simulationParameterLocations))
            for i in range(len(simulationParameterLocations)):
                #print(simulationParameterLocations[i])
                #(gsdFileLocations[i-2])

                self.name = str(gsdFileLocations[i]).split('/')[-1].split('.')[0]
                location = str(gsdFileLocations).split('/')[0] + "/" + str(gsdFileLocations).split('/')[1] + "/" + str(gsdFileLocations).split('/')[2] + "/ "

                self.getSimulationParameters(simulationParameterLocations[i])
                self.constructPolymerObjects(gsdFileLocations[i])
                self.constructGeneralPolymerProfiles(saveLocation = location)
                self.getPositionProbabilityData(name = self.name + "_pos_density",location=location)
                self.animatePositionData(fps=100)
        def animatePositionData(self,fps=500,Interval = 0):
            fig, ax = plt.subplots()
            fig.set_tight_layout(True)
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
            artists = []
            for i in range(len(self.gsd_data[0].particles.position)):
                drawer = plt.Circle([0,0],radius=0.3,color='b',linewidth=0.5,fill=False, clip_on = False)
                artists.append(drawer)

            def update(i):

                # Update the line and the axes (with a new xlabel). Return a tuple of
                # "artists" that have to be redrawn for this frame.

                x = self.gsd_data[int(i)].particles.position[:,0]
                y = self.gsd_data[int(i)].particles.position[:,1]
                ax.clear()
                if circles == True:
                    for j in range(len(x)):
                        artists[j].center = (x[j],y[j])
                        ax.add_patch(artists[j])
                for i in range(self.parameters.getNumberChains()):
                    length = self.parameters.getLength()
                    ax.plot(x[i*length:i*length + length],y[i*length:i*length + length])
                    ax.plot(x[i*length:i*length + length],y[i*length:i*length + length],'o')

                ax.set_xlim(xmin,xmax)
                ax.set_ylim(ymin,ymax)
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

            #print(self.parameters.getNumberChains(),self.parameters.getLength(),len(self.gsd_data))
            for i in range(int(self.parameters.getNumberChains())):

                x = []
                y = []
                for j in range(int(self.parameters.getLength())):

                    x.append([])
                    y.append([])
                    for k in range(int(len(self.gsd_data)*self.interval),len(self.gsd_data)):

                        x[-1].append(self.gsd_data[k].particles.position[i*self.parameters.getNumberChains() + j,0])
                        y[-1].append(self.gsd_data[k].particles.position[i*self.parameters.getNumberChains() + j,1])

                self.polymers.append(PolymerObject([x,y]))


        def gaussian(self,x, mu, sig):
            return numpy.exp(-numpy.power(x - mu, 2.) / (2 * numpy.power(sig, 2.)))

        def constructGeneralPolymerProfiles(self,savePlot = True,saveLocation=""):

            if savePlot == True:
                for i in range(len(self.polymers)):
                    plt.clf()

                    generalProfileWidths = self.polymers[i].getWidths()
                    #y = range(len(generalProfileWidths),0,-1)
                    y = range(len(generalProfileWidths))
                    mu = 0
                    print(mu)
                    #x =  numpy.linspace(-1 + mu, 1 + mu, 120)
                    x = [mu]*len(y)
                    xerr = []
                    for j in range(len(y)):
                        xerr.append(generalProfileWidths[j][0])
                    plt.errorbar(x,y,xerr=xerr)
                    #for j in range(len(y)):
                    #
                        #plt.plot(x, self.gaussian(x, mu , generalProfileWidths[j][0]) + 3*j,color='k',linewidth = 1)
                    plt.savefig(self.name + "_p" + str(i)+"_widths.png")
                    os.system("mv " + self.name + "_p" + str(i)+ ".png " + saveLocation)






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
            plt.savefig(name + ".png")
            os.system("mv " +name+".png " + location )




class PolymerObject():
    def __init__(self,data):
        self.particle = data
        self.generateWidths()
        self.calculateTilt()
    def getParticleData(self,id):
        return [self.particle[0][id],self.particle[1][id]]

    def getParticleWidth(self,id):
        return self.particleWidth[id]

    def generateWidths(self):
        self.particleWidth = []
        for i in range(len(self.particle[0])):
            self.particleWidth.append([   self.stdDev(self.particle[0][i])  , self.stdDev(self.particle[1][i])])

    def calculateTilt(self):
        self.tilt = 0
        self.mean = sum(self.particle[0][0])/len(self.particle[0][0])
        rise = 0
        run =0
        for i in range(len(self.particle[0][0])):
            run += self.particle[0][-1][i]
            run -= self.particle[0][0][i]
            rise += self.particle[1][-1][i]
            rise -= self.particle[1][0][i]

        self.tilt = run/rise * 1/len(self.particle[0][0])

    def getTilt(self):
        return self.tilt
    def getMean(self):
        return self.mean
    def getWidths(self):
        return self.particleWidth

    def stdDev(self,data):
        std = 0
        mean = sum(data)/len(data)
        difference = 0
        for i in range(len(data)):
            difference += (data[i]-mean)**2
        std = math.sqrt(1/(len(data)-1)*difference)

        return std
