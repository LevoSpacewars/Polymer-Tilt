
import os
import hoomd
import hoomd.md
import math
import random
import numpy
import gsd
import gsd.hoomd
def hardContact(r, rmin, ramx, l_0,K):
    V = K*(r-l_0)**2
    F = -2*K*(r-l_0)
    return (V,F)

def harmonicF(r, rmin, rmax,l_0,l_max,K):
    if r < l_0:
        V = K*(r-l_0)**2
        F = -2*K*(r-l_0)
        return (V,F)

    V = - K*l_max**2/2*math.log(1-(r-l_0)**2/l_max**2)
    F = - (K*(r-l_0))/(1-(r-l_0)**2/l_max**2)
    return (V,F)


class PolymerSimulationParameters():
    def __init__(self,data = None):
        if data == None:
            print("********\n*\n*NO PARAMETERS HAVE BEEN PASSED INTO POLYMER_SIMULATION_PARAMETERS\n*initializing all parameters to 0\n*\n********")
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
            a = 0
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
    def __init__(self,name=None,parameter=None,initializer='--mode=cpu'):
        hoomd.context.initialize(initializer)
        self.parameter = parameter
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

    def run(self):
        self.setupFileSystem()
        self.simulationReadMeDump()
        hoomd.analyze.log(filename="energy.log",
                          quantities=['potential_energy', 'temperature'],
                          period=self.parameter.getProbePeriod(),
                          overwrite=True);
        gsdname="polymer.gsd"
        hoomd.dump.gsd(gsdname, period=self.parameter.getProbePeriod(), group=self.all, overwrite=True);



        hoomd.run(self.parameter.getRunLength())
        os.system("mv polymer.gsd " + self.DirectoryName + "/" + self.currentSimulation)
        os.system("mv energy.log " +self.DirectoryName + "/" + self.currentSimulation)
        return self.DirectoryName + "/" + self.currentSimulation

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
        lines       = self.parameter.getNumberChains()
        length      = self.parameter.getLength()
        l_0         = self.parameter.getPairRadiusEqualibrium()
        l_max       = self.parameter.getPairMaximumRadius()
        K           = self.parameter.getPairPotentialStrength()
        amplitude   = self.parameter.getPeriodicAmplitude()
        pull        = self.parameter.getPullForce()
        rmax        = l_0 + l_max - 0.000001

        self.all    = hoomd.group.all()
        self.anchor = hoomd.group.type(name='anchor', type='A')
        self.pulley = hoomd.group.type(name='pulley',type='B')
        self.chain  = hoomd.group.type(name='chain', type='C')

        nl = hoomd.md.nlist.cell();



        table = hoomd.md.pair.table(width=100000, nlist=nl) #I think the radius cut is the ran
        table.pair_coeff.set('A', 'A', func=hardContact,rmin=0,rmax=l_0, coeff=dict(l_0=l_0, K=K))
        table.pair_coeff.set('A', 'B', func=hardContact,rmin=0,rmax=l_0, coeff=dict(l_0=l_0, K=K))
        table.pair_coeff.set('B', 'B', func=hardContact,rmin=0,rmax=l_0, coeff=dict(l_0=l_0, K=K))

        table.pair_coeff.set('B', 'C', func=hardContact,rmin=0,rmax=l_0, coeff=dict(l_0=l_0, K=K))
        table.pair_coeff.set('A', 'C', func=hardContact,rmin=0,rmax=l_0, coeff=dict(l_0=l_0, K=K))
        table.pair_coeff.set('C', 'C', func=hardContact,rmin=0,rmax=l_0, coeff=dict(l_0=l_0, K=K))

        nl.reset_exclusions(exclusions = [])

        harmonic = hoomd.md.bond.table(width = 10000000);
        harmonic.bond_coeff.set('polymer', func = harmonicF,rmin=0, rmax=rmax,coeff = dict(l_0=l_0,l_max=l_max,K=K));



        self.tensionForce = hoomd.md.force.constant(group = self.pulley, fvec=(0.0,pull,0.0))  # ?
        periodic = hoomd.md.external.periodic()
        periodic.force_coeff.set('A', A=amplitude, i=0, w=1, p=lines)
        periodic.force_coeff.set('B', A=amplitude, i=0, w=1, p=lines)
        periodic.force_coeff.set('C', A=amplitude, i=0, w=1, p=lines)

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


    def simulationReadMeDump(self):
        text = open(self.DirectoryName + "/" + self.currentSimulation + "/simulation_parameters.txt","w+")
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
        def __init__(self,directory):
            self.directory = directory

            self.getSimulationParameters(directory)
            self.constructPolymerObjects()
            self.constructGeneralPolymerProfile()

        def getSimulationParameters(self,directory):
            file = open(directory + "/simulation_parameters.txt",'r');
            data = []
            lines = file.readlines()
            for i in range(len(lines)-1):

                data.append(float(lines[i].split('=')[1]))
            data.append([float(lines[len(lines)-1].split('=')[1].split(',')[0]),float(lines[len(lines)-1].split('=')[1].split(',')[1])])

            self.parameters = PolymerSimulationParameters(data=data)


        def constructPolymerObjects(self):

            self.gsd_data =gsd.hoomd.open(self.directory  + "/polymer.gsd",'rb')
            self.polymers = []

            print(self.parameters.getNumberChains(),self.parameters.getLength(),len(self.gsd_data))
            for i in range(int(self.parameters.getNumberChains())):

                x = []
                y = []
                for j in range(int(self.parameters.getLength())):

                    x.append([])
                    y.append([])
                    for k in range(len(self.gsd_data)):

                        x[-1].append(self.gsd_data[k].particles.position[i*self.parameters.getNumberChains() + j,0])
                        y[-1].append(self.gsd_data[k].particles.position[i*self.parameters.getNumberChains() + j,1])

                self.polymers.append(PolymerObject([x,y]))

        def constructGeneralPolymerProfile(self,savePlot = True):
            self.generalProfileWidths = [[],[]]
            for i in range(self.parameters.getLength()):
                self.generalProfileWidths[0].append(0)
                self.generalProfileWidths[1].append(0)

            x=[]
            for i in range(len(self.polymers)):
                for j in range(self.parameters.getLength()):
                    self.generalProfileWidths[0][j] += (self.polymers[i].getParticleWidth(j)[0])**2
                    self.generalProfileWidths[1][j] += (self.polymers[i].getParticleWidth(j)[1])**2



            for j in range(self.parameters.getLength()):
                self.generalProfileWidths[0][j] = math.sqrt(self.generalProfileWidths[0][j])
                self.generalProfileWidths[1][j] = math.sqrt(self.generalProfileWidths[1][j])
            print(self.generalProfileWidths)
            if savePlot == True:
                from matplotlib import pyplot as plt

                y = range(0,len(self.generalProfileWidths[0])*2,2)
                x = [0]*(len(y))
                xerr = []
                yerr = []
                for i in range(len(self.generalProfileWidths[0])):
                    xerr.append(self.generalProfileWidths[0][i])
                    yerr.append(self.generalProfileWidths[1][i])
                print(len(x),len(y),len(xerr),len(yerr))
                plt.errorbar(x,y,xerr=xerr,yerr=yerr)
                plt.savefig("test.png")
class PolymerObject():
    def __init__(self,data):
        self.particle = data
        self.generateWidths()
    def getParticleData(self,id):
        return [self.particle[0][id],self.particle[1][id]]

    def getParticleWidth(self,id):
        return self.particleWidth[id]

    def generateWidths(self):
        self.particleWidth = []
        for i in range(len(self.particle[0])):
            self.particleWidth.append([   self.stdDev(self.particle[0][i])  , self.stdDev(self.particle[1][i])])
        from matplotlib import pyplot as plt
        x= []
        y=[]
        for i in range(len(self.particle[0])):
            x.append(self.particle[0][i][-1])
            y.append(self.particle[1][i][-1])
        plt.plot(x,y,'.')
        plt.savefig("pos.png")


    def stdDev(self,data):
        std = 0
        mean = sum(data)/len(data)
        difference = 0
        for i in range(len(data)):
            difference += (data[i]-mean)**2
        std = math.sqrt(1/(len(data)-1)*difference)

        return std
