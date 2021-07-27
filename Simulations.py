
import os
from typing import Tuple, List
import hoomd
import hoomd.md
import math
import random
import numpy
import gsd
import gsd.hoomd
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation

#from PIL import Image
circles = False


def hardContact(r, rmin, ramx, l_0,K):
    V = 1/2 *K*(r-l_0)**2
    F = -K*(r-l_0)
    return (V,F)

def particlePotential(r, rmin, rmax,paricle_radius,max_bond_radius,strength_coef):
    # r                     : float :   radius between two particles
    # rmin                  : float :   minimum radius
    # rmax                  : float :   maximum radius
    # paricle_radius        : float :   The radius of each particle (also the equlibrium point)
    # max_bond_radius       : float :   the maximum distance between the surfaces of 2 particles
    # strength_coef         : float :   Prefactor meant to increase the effectiveness of the function


    contactRegime = paricle_radius + 0.00001
    correctionRegime = paricle_radius + max_bond_radius
    particleDistance = r - paricle_radius
    correction = max_bond_radius - 0.0001

    if r < contactRegime:
        potential = strength_coef * (particleDistance)**2
        force = -2 * strength_coef * (particleDistance)
        return (potential,force)


    elif r >= correctionRegime:
        potential = - 0.5 * strength_coef * max_bond_radius**2 * math.log( 1 - correction**2 / max_bond_radius**2 )
        force = - ( strength_coef * correction ) / ( 1 - correction**2 / max_bond_radius**2 )
        return (potential,force)


    else:
        potential = - 0.5 * strength_coef * max_bond_radius**2 * math.log( 1- particleDistance**2 / max_bond_radius**2 )
        force =  - ( strength_coef * particleDistance ) / ( 1 - particleDistance**2 / max_bond_radius**2 )
        return (potential,force)
def particlePotentialHarmonic(r, rmin, rmax,particle_diameter,strength_coef):
    # r                     : float :   radius between two particles
    # rmin                  : float :   minimum radius
    # rmax                  : float :   maximum radius
    # particle_diameter     : float :   The radius of each particle (also the equlibrium point)
    # strength_coef         : float :   Prefactor meant to increase the effectiveness of the function


    particleDistance = r - particle_diameter

    potential = 1/2*(strength_coef * particleDistance**2)
    force     = -strength_coef * particleDistance
    return (potential, force)

class DisorderParameter():
    def __init__(self, random_seed, amplitude_range,nodes, width) -> None:
        self.bond_radius = 0.05
        self.width = width
        self.seed = random_seed
        self.amplitude_range = amplitude_range
        self.nodes = nodes
        self.disorder:List[Tuple[float, int, float]]  = self.generate_disorder()

    def generate_disorder(self)-> List[Tuple[float, int, float]] :
        import random as r

        disorder_level = int(self.nodes)
        seed = self.seed
        r.seed(seed)
        inv_lamda = int(1/(self.bond_radius))
        Arange = self.amplitude_range

        p:List[Tuple[float, int, float]] = []  # (amplitude, nodes, phase)

        for null in range(disorder_level):
            A = r.random() * (Arange[1] - Arange[0]) + Arange[0]
            phase = 2 * r.random() * math.pi
            nodes = r.randint(1, inv_lamda)
            p.append((A, nodes, phase))
        return p

    def get_disorder(self):
        return self.disorder

    def write_disorder(self,path):
        wfile = open(path + '/disorder_info.txt','w')
        wfile.write("Amplitude, nodes, Phase\n")
        for i in self.disorder:
            wfile.write(f"{i[0]},{i[1]},{i[2]}\n")
        wfile.close()

    def import_disorder(self,path):
        import numpy as np
        array = List(np.genfromtxt(path + "/disorder_info.txt", delimiter=',',skip_header=1))
        return array

class PolymerSimulationParameters():
    def __init__(self,sheerforcerange=[0,0],df=0,length=0,lines=0,rez=0,K=0,l_0=0,l_max=0,pull=0,amplitude=0,gamma=0,kbT=0,dt=0,probePeriod=0,runLength=0,boxdimx=0,boxdimy=0,integraor = "brownian", direction = "forward", disorder = None):
        self.sheerForceRange=sheerforcerange
        self.df=df
        self.length=int(length)
        self.lines=int(lines)
        self.rez=rez
        self.K=K
        self.l_0=l_0
        self.l_max=l_max
        self.pull=pull
        self.amplitude=amplitude
        self.gamma=gamma
        self.kbT=kbT
        self.dt=dt
        self.probePeriod=int(probePeriod)
        self.runLength=int(runLength)
        self.boxdimx= boxdimx
        self.boxdimy=boxdimy
        self.integrator="brownian"
        self.direction = direction

    def setSheerForceRange(self,x1,x2):
        self.sheerForceRange = [x1,x2]
    def setDf(self,x):
        self.df = x
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
    def setBoxDimx(self,x):
        self.boxdimx = x
    def setBoxDimy(self,x):
        self.boxdimy = x
    def setIntegrator(self,x):
        self.integrator = x
    def setRunDirection(self,x):
        self.direction  =x

    def getSheerForceRange(self):
        return self.sheerForceRange
    def getDf(self):
        return self.df
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
    def getBoxDimx(self):
        return self.boxdimx
    def getBoxDimy(self):
        return self.boxdimy
    def getIntegrator(self):
        return self.integrator
    def getRunDirection(self):
        return self.direction

    def writeParameters(self,name="",dir=""):
        text = None
        text = open(dir +"/" + name + "_simulation_parameters.txt","w+")

        text.write("sheerForceRange=" + str(self.getSheerForceRange()[0]) + "," + str(self.getSheerForceRange()[1])       + "\n")
        text.write("df=" + str(self.getDf())                                 + "\n")
        text.write("length="    +      str(self.getLength())                 + "\n")
        text.write("lines="     +      str(self.getNumberChains())           + "\n")
        text.write("rez="       +      str(self.getPairRadius())             + "\n")
        text.write("K="         +      str(self.getPairPotentialStrength())  + "\n")
        text.write("l_0="       +      str(self.getPairRadiusEqualibrium())  + "\n")
        text.write("l_max="     +      str(self.getPairMaximumRadius())      + "\n")
        text.write("pull="      +      str(self.getPullForce())              + "\n")
        text.write("amplitude=" +      str(self.getPeriodicAmplitude())      + "\n")
        text.write("gamma="     +      str(self.getGamma())                  + "\n")
        text.write("kbT="       +      str(self.getKBT())                    + "\n")
        text.write("dt="        +      str(self.getTimeStep())               + "\n")
        text.write("probePeriod="+     str(self.getProbePeriod())            + "\n")
        text.write("runLength=" +      str(self.getRunLength())              + "\n")
        text.write("BoxDimx="    +      str(self.getBoxDimx())               + "\n")
        text.write("BoxDimy="    +      str(self.getBoxDimy())               + "\n")
        text.write("Integrator=" +     str(self.getIntegrator())             + "\n")
        text.write("Direction=" +       str(self.getRunDirection())          + "\n")
        text.close()
    def loadParameters(self,fileLocation):
        file = open(fileLocation,'r');
        lines = file.readlines()
        for i in range(len(lines)):
            obj = lines[i]
            if "sheerForceRange=" in obj:
                self.setSheerForceRange(float(lines[i].split('=')[1].split(',')[0]),float(lines[i].split('=')[1].split(',')[1]))
            elif "df=" in obj:
                self.setDf(int(lines[i].split('=')[1]))
            elif "length=" in obj:
                self.setLength(int(lines[i].split('=')[1]))
            elif "lines=" in obj:
                self.setNumberChains(int(lines[i].split('=')[1]))
            elif "rez=" in obj:
                self.setPairRadius(float(lines[i].split('=')[1]))
            elif "l_0=" in obj:
                self.setPairRadiusEqualibrium(float(lines[i].split('=')[1]))
            elif "l_max=" in obj:
                self.setPairMaximumRadius(float(lines[i].split('=')[1]))
            elif "pull=" in obj:
                self.setForcePull(float(lines[i].split('=')[1]))
            elif "amplitude=" in obj:
                self.setPeriodicAmplitude(float(lines[i].split('=')[1]))
            elif "gamma=" in obj:
                self.setGamma(float(lines[i].split('=')[1]))
            elif "kbT=" in obj:
                self.setKBT(float(lines[i].split('=')[1]))
            elif "dt=" in obj:
                self.setTimeStep(float(lines[i].split('=')[1]))
            elif "probePeriod=" in obj:
                self.setProbePeriod(int(lines[i].split('=')[1]))
            elif "runLength=" in obj:
                self.setRunLength(int(lines[i].split('=')[1]))
            elif "BoxDimx=" in obj:
                self.setBoxDimx(float(lines[i].split('=')[1]))
            elif "BoxDimy=" in obj:
                self.setBoxDimy(float(lines[i].split('=')[1]))
            elif "Integrator=" in obj:
                self.setIntegrator(lines[i].split('=')[1])
            elif "Direction=" in obj:
                self.setRunDirection(lines[i].split('=')[1])
            elif "DisorderLevel" in obj:
                self.setDisorderLevel(lines[i].split('=')[1])
            elif "RandomSeed" in obj:
                self.setRandomSeed(lines[i].split('=')[1])
        file.close()



class PolymerSimulation():
    def __init__(self):
        self.parameter = None
        random.seed()

    def init(self, parameter: PolymerSimulationParameters =None,initializer='--mode=gpu',loadSave=None,dirName=None,probe = False):
        hoomd.context.initialize(initializer)
        self.parameter = parameter
        if loadSave is None:
            random.seed()

            if probe !=False:
                self.setupFileSystem()
            #1. create new directory for this run

            self.populateSystem()
            #2. load and populate snapshot for simulation

            self.initializeForces()
            self.initializeIntegrator()

        else:
            hoomd.init.read_snapshot(hoomd.data.gsd_snapshot(loadSave)) # only used when I am loading a save state
            if probe !=False:
                self.setupFileSystem()
            self.initializeForces()
            self.initializeIntegrator()

    def probe(self, run_id,sheer_value,path,server = False):
        name = str(run_id) + "_sheer_" + str(sheer_value)
        self.setupFileSystem(name=name)
        self.apply_disorder()
        self.view_potential()

        #self.view_potential(self.parameter.getNumberChains())
        nameg = str(run_id) + "_sheer_" + str(sheer_value)+".gsd"
        hoomd.dump.gsd(filename=self.DirectoryName + "/" +"trajectory.gsd", period=self.parameter.getProbePeriod(), group=self.all, overwrite=True)
        self.simulationReadMeDump(singular = True,sheer_value=sheer_value)


        self.tensionForce.set_force(fvec=(sheer_value,self.parameter.getPullForce(),0.0))
        self.sheerForce.set_force(fvec=(-sheer_value,0,0))
        hoomd.run(self.parameter.getRunLength())

        if server == True:
            os.system("mv " + self.DirectoryName + " /projects/softmatter/apatapof/runs/" + self.DirectoryName)

        return path


    def run(self, server = False): #main run function

        self.setupFileSystem()
        print(self.DirectoryName)

        #self.view_potential(self.parameter.getNumberChains())


        file = hoomd.dump.gsd(filename= self.DirectoryName + "/trajectory.gsd", period=self.parameter.getProbePeriod(), group=self.all, overwrite=True)
        hoomd.analyze.log(filename=self.DirectoryName + "/Energy.log",quantities=['potential_energy', 'temperature'],period=self.parameter.getProbePeriod(),overwrite=True);

        self.simulationReadMeDump()

        for i in range(self.parameter.getDf()):



            sheerforce = (i/self.parameter.getDf())*(self.parameter.getSheerForceRange()[1] - self.parameter.getSheerForceRange()[0]) + self.parameter.getSheerForceRange()[0]

            self.tensionForce.set_force(fvec=(sheerforce,self.parameter.getPullForce(),0.0))
            self.sheerForce.set_force(fvec=(-sheerforce,0,0))

            print(sheerforce)

            hoomd.run(self.parameter.getRunLength())

            # creating frame 0 load state
            hoomd.dump.gsd(filename="save_ " + str(sheerforce) + ".gsd", period=None, group=self.all, overwrite=True, dynamic=['attribute', 'property', 'momentum', 'topology'])
            hoomd.run(1)
            os.system("mv " + "save_\ " + str(sheerforce) + ".gsd " + self.DirectoryName + "/")
        #os.system("mv " + "trajectory.gsd" + " " + self.DirectoryName + "/")
        #os.system("mv " + "Energy.log"+ " " +self.DirectoryName + "/")
        os.system("DataHandler/Executable " + " /projects/softmatter/apatapof/jobs/" + self.DirectoryName + "/")
        if server == True:
            os.system("/projects/softmatter/apatapof/Polymer-Tilt/DataHandler/Executable " + " /projects/softmatter/apatapof/jobs/" + self.DirectoryName + "/")
            os.system("mv " + self.DirectoryName + " /projects/softmatter/apatapof/runs/" + self.DirectoryName)





        return self.DirectoryName + "/"


    def revRun(self,server=False): #unused
        self.setupFileSystem()


        hoomd.dump.gsd(filename="trajectory.gsd", period=self.parameter.getProbePeriod(), group=self.all, overwrite=True)
        hoomd.analyze.log(filename="Energy.log",quantities=['potential_energy', 'temperature'],period=self.parameter.getProbePeriod(),overwrite=True);

        self.parameter.setRunDirection("reverse")
        self.simulationReadMeDump()

        for k in range(1,self.parameter.getDf()+1):
            i = self.parameter.getDf() - k
            print(i)
            sheerforce = i/self.parameter.getDf()*(self.parameter.getSheerForceRange()[1] - self.parameter.getSheerForceRange()[0]) + self.parameter.getSheerForceRange()[0]
            print(sheerforce)
            self.tensionForce.set_force(fvec=(sheerforce,self.parameter.getPullForce(),0.0))
            self.sheerForce.set_force(fvec=(-sheerforce,0,0))
            print(sheerforce)
            hoomd.run(self.parameter.getRunLength())



        os.system("mv " + "trajectory.gsd" + " " + self.DirectoryName + "/")
        os.system("mv " + "Energy.log"+ " " +self.DirectoryName + "/")
        if server:
            os.system("mv " + self.DirectoryName + " /projects/softmatter/apatapof/runs")

        return self.DirectoryName + "/"


    def initializeIntegrator(self):

        hoomd.md.integrate.mode_standard(dt=self.parameter.getTimeStep());
        TRamp = []

        for i in range(self.parameter.getDf()): #temperature parameters for annealing in hoomd.variant()
            TRamp.append((int(self.parameter.runLength*i),float(self.parameter.kbT * 0.5 + self.parameter.kbT)))
            TRamp.append( (int(self.parameter.runLength*i + self.parameter.runLength/2), float(self.parameter.kbT)))
            TRamp.append( (int(self.parameter.runLength*i + self.parameter.runLength -1), float(self.parameter.kbT)))
        for i in range(len(TRamp)):
            print(str(type(TRamp[i][0])) + "," + str(type(TRamp[i][1])))


        if (self.parameter.getIntegrator() == "brownian"): # not used, since I've only been using the legavin integrator
            if(self.parameter.getNumberChains() == 1):
                self.bs = hoomd.md.integrate.brownian(group=self.most, kT= hoomd.variant.linear_interp(points = TRamp), seed=random.randint(0,999999));

            else:
                self.bs = hoomd.md.integrate.brownian(group=self.all, kT= hoomd.variant.linear_interp(points = TRamp), seed=random.randint(0,999999));

            self.bs.set_gamma('A', gamma=self.parameter.getGamma())
            self.bs.set_gamma('B', gamma=self.parameter.getGamma())
            self.bs.set_gamma('C', gamma=self.parameter.getGamma())
        else:
            if(self.parameter.getNumberChains() == 1): # for a single polymer all the particles except the base are simulated.
                self.bs = hoomd.md.integrate.langevin(group = self.all, kT= hoomd.variant.linear_interp(points = TRamp), seed=random.randint(0,99999),noiseless_r=True);
            else: #defualt
                self.bs = hoomd.md.integrate.langevin(group = self.all, kT= hoomd.variant.linear_interp(points = TRamp), seed=random.randint(0,99999),noiseless_r=True);

            self.bs.set_gamma('A', gamma=self.parameter.getGamma()) #damping
            self.bs.set_gamma('B', gamma=self.parameter.getGamma())
            self.bs.set_gamma('C', gamma=self.parameter.getGamma())


    def setupFileSystem(self,name=""):

        #  constructs a directory file based on the computer time, for the simulation being run.

        from datetime import datetime
        time = datetime.now().strftime('%H-%M-%S')

        if name != "":
            self.DirectoryName = name
            os.system("mkdir " + name)

        elif self.parameter.getNumberChains() is not 1:
            length = self.parameter.getLength();
            T      = self.parameter.getKBT();
            A      = self.parameter.getPeriodicAmplitude();

            self.DirectoryName = str(length) + "," + str(T) + "," + str(A) + "," + time;
            os.system("mkdir " + self.DirectoryName)
        else:
            length = self.parameter.getLength();
            T      = self.parameter.getKBT();
            A      = self.parameter.getPeriodicAmplitude();
            self.DirectoryName = "BASE," + str(length) + "," + str(T) + "," + str(A) + "," + time;
            os.system("mkdir " + self.DirectoryName)


    def initializeForces(self): #where particle-particle and bond interactions are defined
        added = 0
        lines       = self.parameter.getNumberChains()
        length      = self.parameter.getLength()
        bond_length = self.parameter.getPairRadiusEqualibrium()*2
        K           = self.parameter.getPairPotentialStrength()
        amplitude   = self.parameter.getPeriodicAmplitude()
        pull        = self.parameter.getPullForce()
        width       = self.parameter.getBoxDimx()

        #defining group names for particle types
        self.all    = hoomd.group.all()
        self.anchor = hoomd.group.type(name='anchor', type='A') #anchor is the base particle for each polymer
        self.pulley = hoomd.group.type(name='pulley',type='B') # pulley is the top particle for each polymer
        self.chain  = hoomd.group.type(name='chain', type='C') # everything else
        self.most   = hoomd.group.union(name='most',  a = self.chain, b = self.pulley) #everthing but the anchor group

        nl = hoomd.md.nlist.cell();



        table = hoomd.md.pair.table(width=10000, nlist=nl) # describes stiff quadratic potential for particle-particle interactions
        table.pair_coeff.set('A', 'A', func=hardContact,rmin=0,rmax=bond_length, coeff=dict(l_0=bond_length, K=K))
        table.pair_coeff.set('A', 'B', func=hardContact,rmin=0,rmax=bond_length, coeff=dict(l_0=bond_length, K=K))
        table.pair_coeff.set('B', 'B', func=hardContact,rmin=0,rmax=bond_length, coeff=dict(l_0=bond_length, K=K))

        table.pair_coeff.set('B', 'C', func=hardContact,rmin=0,rmax=bond_length, coeff=dict(l_0=bond_length, K=K))
        table.pair_coeff.set('A', 'C', func=hardContact,rmin=0,rmax=bond_length, coeff=dict(l_0=bond_length, K=K))
        table.pair_coeff.set('C', 'C', func=hardContact,rmin=0,rmax=bond_length, coeff=dict(l_0=bond_length, K=K))

        nl.reset_exclusions(exclusions = [])

        harmonic = hoomd.md.bond.table(width = 100000); #defining the bond potential
        harmonic.bond_coeff.set('polymer', func = particlePotentialHarmonic,rmin=0, rmax=100,coeff = dict(particle_diameter=bond_length, strength_coef=K)); #bond potential


        self.tensionForce = hoomd.md.force.constant(group = self.pulley, fvec=(0.0,0,0.0))  # FORCES INTIALIZED HERE -> CHANGED IN RUN FUNCTION
        self.sheerForce   = hoomd.md.force.constant(group = self.anchor, fvec=(0.0,0.0,0.0))
        periodic = hoomd.md.external.periodic() #External potential defined
        if lines is not 1:
            periodic.force_coeff.set('A', A=amplitude, i=0, w=1, p=lines+added)
            periodic.force_coeff.set('B', A=amplitude, i=0, w=1, p=lines+added)
            periodic.force_coeff.set('C', A=amplitude, i=0, w=1, p=lines+added)
            print("multipolymer settings")
        else:
            periodic.force_coeff.set('A', A=amplitude, i=0, w=1, p=width)
            periodic.force_coeff.set('B', A=amplitude, i=0, w=1, p=width)
            periodic.force_coeff.set('C', A=amplitude, i=0, w=1, p=width)
            print("single polymer ssetting")
            #periodic.force_coeff.set('A', A=-10000000.0, i=0, w=0, p=10)


        periodic.force_coeff.set('A', A=10000000.0, i=1, w=0, p=500) #used to keep anchor on y=0

    def set_disorder(self, random_seed, amplitude_range,nodes, width):
        self.disorder = DisorderParameter(random_seed, amplitude_range,nodes,width)

    def apply_disorder(self):
        dparam = self.disorder.get_disorder()
        print(len(dparam))
        for tup in dparam:

            f_disorder = hoomd.md.external.periodic()
            f_disorder.force_coeff.set('A', A=tup[0], i=0, w=tup[2], p=tup[1])
            f_disorder.force_coeff.set('B', A=tup[0], i=0, w=tup[2], p=tup[1])
            f_disorder.force_coeff.set('C', A=tup[0], i=0, w=tup[2], p=tup[1])

        self.disorder.write_disorder(self.DirectoryName)


    def view_potential(self):
        print("generating potential")
        import math as m
        width = self.parameter.lines
        cos_param = self.disorder.get_disorder()
        #print(cos_param)
        #cos_param = []
        cos_param.append((self.parameter.amplitude,width, width, 0))
        N = 1000
        import numpy as np
        from matplotlib import pyplot as plt


        den = width / N
        potential = np.zeros(N * N).reshape(N,N)
        for element in cos_param:
            for x in range(len(potential[0])):
                p = x * den - width/2
                A = element[0]
                freq = element[1] / width * 2 * math.pi
                phase = element[2]
                potential[0][x] += A * m.cos(p * freq  + phase)


        sm = 0
        for i in range(len(cos_param)-1):
            element = cos_param[i]
            sm += element[0] * element[0]
        sm = math.sqrt(sm)
        fm = abs(sm/ self.parameter.amplitude)
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)


        from matplotlib import pyplot as plt
        fig, axs = plt.subplots(2)
        axs[0].imshow(potential)
        x = np.linspace(-width/2,width/2,len(potential[0]))
        axs[1].plot(x,potential[0])
        plt.text(-width/2,max(potential[0]),f"disorderAmplitude/periodicAmplitude= {round(fm,4)}",bbox=props,fontsize=6,verticalalignment='top');
        plt.savefig(self.DirectoryName + "/potential.png", dpi=300, bbox_inches='tight')
        plt.close()


    def definePositions(self):
        lines   = self.parameter.getNumberChains()
        length  = self.parameter.getLength()
        l_0     = self.parameter.getPairRadiusEqualibrium()
        l_max   = self.parameter.getPairMaximumRadius()
        pos     = []

        if lines == 1:
            x = 0
            y = 0
            pos.append([x,y,0])
            for j in range(1,length):
                y = y + 2*l_0
                pos.append([x,y,0])

        else:

            for i in range(lines):
                x = i - lines/2
                y = 0
                pos.append([x,y,0])
                for j in range(1,length):
                    y = y + 2*l_0
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

    def defineDiameter(self):
        d = []
        lines   = self.parameter.getNumberChains()
        length  = self.parameter.getLength()
        for i in range(lines*length):
            d.append(2*self.parameter.getPairRadiusEqualibrium())
        return d

    def populateSystem(self):
        pos = []
        diameter = []
        bonds = []
        id = []
        mass = []
        types = ['A','B','C']

        self.boxdim=[0,0]
        self.boxdim[0] = self.parameter.getNumberChains()
        self.boxdim[1] = 1000
        self.parameter.setBoxDimx(self.boxdim[0])
        self.parameter.setBoxDimy(self.boxdim[1])


        if self.parameter.getNumberChains() == 1:
            self.boxdim[0] = 10;
            self.parameter.setBoxDimx(self.boxdim[0])


        boundary = hoomd.data.boxdim(Lx = self.boxdim[0], Ly = self.boxdim[1], dimensions=2)
        snapshot = hoomd.data.make_snapshot(N=  int(self.parameter.getNumberChains()  *  self.parameter.getLength()), box=boundary, particle_types=types, bond_types=['polymer'])

        pos     = self.definePositions()
        bonds   = self.defineBonds()
        id      = self.defineIds()
        mass    = self.defineMass()
        diameter= self.defineDiameter()

        snapshot.particles.position[:] = pos
        snapshot.bonds.resize(len(bonds))
        snapshot.bonds.group[:] = bonds
        snapshot.particles.typeid[:] = id
        snapshot.particles.mass[:] = mass
        snapshot.particles.diameter[:] = diameter

        hoomd.init.read_snapshot(snapshot)

    def simulationReadMeDump(self,name="", dir="", singular = False, sheer_value = 0):
        if singular == True:
            self.parameter.setSheerForceRange(sheer_value,sheer_value + 1)
            self.parameter.setDf(1)
        self.parameter.writeParameters(name=name,dir=self.DirectoryName)

#################################### only for extracting data ####################################################
class DataVisualizer():
        def __init__(self,basedirectory="",interval=0):

            #basedirectory  : string : base directory for the simulation in question
            #interval       : float  : 0.0-1.0 which is used to offset the start of sampling

            self.interval = interval
            self.directory = basedirectory


        def init(self, plotTilt= False, plotPolymerProfiles = False, plotProbabilityMap = False, plotEnergy = False, animatePolymers = False):
            import glob

            location = self.directory

            self.parameters = None
            self.polymers = []
            self.forceValues = []



            self.name = "trajectory.gsd"
            location = location.split('/')[0] + "/"
            self.gsd_data =gsd.hoomd.open(location + self.name,'rb')
            print(self.gsd_data[0])
            print(len(self.gsd_data))
            parameterfilename=location+ "_simulation_parameters.txt"
            dataname= location + self.name

            self.getSimulationParameters(parameterfilename)
            self.forceValues = self.getForceRange()
            if plotProbabilityMap:
                self.plotPositionProbabilityData(name = self.name,location=dataname)
            if plotTilt or plotPolymerProfiles:
                self.constructPolymerObjects(dataname)

            if plotTilt:
                self.plotTilt()
                self.plotLengthDx()
            if plotPolymerProfiles:
                self.plotGeneralPolymerProfiles()

            if plotEnergy:
                self.plotEnergy(location="ADDNAME",name=self.name+"energy.png")
            if animatePolymers:
                self.animatePositionData(dataname,fps=200,Interval = self.interval)

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
                    p[-1].append(math.sin(j*2*math.pi/rez*boxdimx - boxdimx/2))

            return p



        def animatePositionData(self,fileLocation,fps=250,Interval = 0):

            circles = False
            print("begining animation")
            fig, ax = plt.subplots(figsize=(10, 10))
            # xmax = self.gsd_data[0].particles.position[:,0][0]
            # xmin = xmax
            # ymax = self.gsd_data[0].particles.position[:,1][0]
            # ymin = ymax
            # for i in range(len(self.gsd_data)):
            #     for j in range(len(self.gsd_data[int(i)].particles.position[:,0])):
                    # if self.gsd_data[int(i)].particles.position[:,0][j] < xmin:
                    #     xmin = self.gsd_data[int(i)].particles.position[:,0][j]
                    # if self.gsd_data[int(i)].particles.position[:,0][j] > xmax:
                    #     xmax = self.gsd_data[int(i)].particles.position[:,0][j]

                    # if self.gsd_data[int(i)].particles.position[:,1][j] < ymin:
                    #     ymin = self.gsd_data[int(i)].particles.position[:,1][j]
                    # if self.gsd_data[int(i)].particles.position[:,1][j] > ymax:
                    #     ymax = self.gsd_data[int(i)].particles.position[:,1][j]
            print("finished sorting. sorting animation")
            potential = self.generatePotential(self.parameters.getNumberChains(),1000)
            def update(i):

                # Update the line and the axes (with a new xlabel). Return a tuple of
                # "artists" that have to be redrawn for this frame.
                x = self.gsd_data[int(i)].particles.position[:,0]
                y = self.gsd_data[int(i)].particles.position[:,1]
                print("timeStep:" + str(i) + "/" + str(len(self.gsd_data)) )

                ax.clear()
                #ax.set_size_inches(10,10)
                ax.imshow(potential,extent=[-self.parameters.getNumberChains()/2,self.parameters.getNumberChains()/2,0,self.parameters.getLength()*0.7])
                for i in range(self.parameters.getNumberChains()):
                    length = self.parameters.getLength()
                    #ax.plot(x[i*length:i*length + length],y[i*length:i*length + length])
                    ax.plot(x[i*length:i*length + length],y[i*length:i*length + length],'o')

                # ax.set_xlim(xmin,xmax)
                # ax.set_ylim(ymin,ymax)
                ax.set_aspect('auto')
                ax.set_xlabel(str(i) + "/" + str(len(self.gsd_data)))

                return ax
            indexrange = int(self.parameters.getRunLength() / self.parameters.getProbePeriod())
            for m in range(len(self.forceValues)):
                print(m)
                anim = FuncAnimation(fig, update, frames=numpy.arange((m+Interval)*indexrange , (m+1)*indexrange), interval=int(fps))
                anim.save(self.name + "_F_" + str(self.forceValues[m]) + "_animation"+'.gif', dpi=80, writer='imagemagick')
            print("finished")


        def getSimulationParameters(self,fileLocation):
            self.parameters = PolymerSimulationParameters()
            print(type(self.parameters));
            self.parameters.loadParameters(fileLocation=fileLocation)

        def getForceRange(self):
            forcerange=[]
            for i in range(self.parameters.getDf()):
                forcerange.append(i/self.parameters.getDf() * (self.parameters.getSheerForceRange()[1] - self.parameters.getSheerForceRange()[0]) + self.parameters.getSheerForceRange()[0])
            print(forcerange)
            return forcerange

        def constructPolymerObjects(self,fileLocation):
            import time
            print("constructing polymer Objects")

            self.forceValues = self.getForceRange()
            timer = 0
            indexrange = int(self.parameters.getRunLength() / self.parameters.getProbePeriod())
            for m in range(len(self.forceValues)):
                self.polymers.append([])
                for i in range(int(self.parameters.getNumberChains())):
                    timer = time.perf_counter()
                    print(str(m) + "." + str(i))
                    x = []
                    y = []
                    for k in range(int(indexrange*(m + self.interval)),indexrange*(m+1)):
                        x.append([])
                        y.append([])
                        for j in range(int(self.parameters.getLength())):
                            x[-1].append(self.gsd_data[k].particles.position[i*self.parameters.getLength() + j,0])

                            y[-1].append(self.gsd_data[k].particles.position[i*self.parameters.getLength() + j,1])


                    self.polymers[-1].append(PolymerObject([x,y],L=self.parameters.getNumberChains()))

                    print('fin');



                    print("polymer construction time: " + str(time.perf_counter() - timer))


            print("done Constructing Polymers")
        def gaussian(self,x, mu, sig):
            return numpy.exp(-numpy.power(x - mu, 2.) / (2 * numpy.power(sig, 2.)))

        def plotGeneralPolymerProfiles(self,savePlot = True):
            print("plotting general profiles")
            if savePlot == True:
                for i in range(len(self.polymers)):
                    for j in range(len(self.polymers[i])):
                        plt.clf()
                        generalProfileWidths = self.polymers[i][j].getWidths()
                        #y = range(len(generalProfileWidths),0,-1)
                        y=[]
                        x=[]
                        mu = self.polymers[i][j].getMean()
                        x = mu[0]
                        y = mu[1]

                        xerr = generalProfileWidths[0]
                        plt.plot(x,y)
                        plt.errorbar(x,y,xerr=xerr,fmt='o')
                        t = numpy.linspace(start=x[0]-2, stop=x[-1]+2, num=200)
                        for k in range(len(x)):
                            plt.plot(t, self.gaussian(t, x[k] , generalProfileWidths[k]) + y[k],color='k',linewidth = 1)
                        print(self.name + "_F_" + str(self.forceValues[i]) + "_N_" + str(j) +"_widths.png")
                        plt.title("Polymer Profile")
                        plt.savefig(self.name + "_F_" + str(self.forceValues[i]) + "_N_" + str(j) +"_widths.png")
                        #os.system("mv " + self.name + "_p" + str(i)+ ".png " + saveLocation)


        def grabParticleFrame(self, interval):
            step = self.parameters.getRunLength()/self.parameters.getProbePeriod()
            for i in range(self.parameters.getDf()):
                location = int(step * interval + step*i)
                print(len(self.gsd_data))
                print(location);
                x = self.gsd_data[location].particles.position[:,0]
                y = self.gsd_data[location].particles.position[:,1]
                p = self.generatePotential(10,1000)
                plt.imshow(p,extent=[-self.parameters.getNumberChains()/2,self.parameters.getNumberChains()/2,0,self.parameters.getLength()*0.2])
                plt.plot(x,y,'.',color='r')
                plt.title("Polymer System")
                plt.xlabel("$\hat{x}$")
                plt.ylabel("$\hat{y}$")
                plt.axes().set_aspect('auto')
                plt.show()

        def plotLengthDx(self):
            self.forceValues = self.getForceRange()
            timer = 0

            length=[]
            ulength=[]

            dx=[]
            udx=[]

            f = []
            output = []
            uoutput = []
            for m in range(len(self.forceValues)):
                f.append(self.forceValues[m]/self.parameters.getPullForce())
                dx.append(0)
                length.append(0)
                udx.append(0)
                ulength.append(0)
                output.append(0)
                uoutput.append(0)
                for j in range(len(self.polymers[m])):
                    dx[-1] += self.polymers[m][j].getDx()
                    length[-1] += self.polymers[m][j].getLength()
                dx[-1] = dx[-1]/len(self.polymers[m])
                length[-1] = length[-1]/len(self.polymers[m])
                output[-1] = dx[-1]/length[-1]
                for j in range(len(self.polymers[m])):
                    udx[-1] += (self.polymers[m][j].getDx() - dx[-1])**2
                    ulength[-1] += (self.polymers[m][j].getLength() - length[-1])**2
                udx[-1] = math.sqrt( udx[-1] * 1/(len(self.polymers[m])-1) )
                ulength[-1] = math.sqrt( ulength[-1] * 1/(len(self.polymers[m])) )
                square = ( udx[-1] / dx[-1] )**2 + ( ulength[-1] / length[-1] )**2
                uoutput[-1] = output[-1] * math.sqrt(square)

            plt.clf()

            plt.title('dX/Length vs SheefForce/Tension: Ft:' + str(self.parameters.getPullForce()) + ' kbT=' + str(self.parameters.getKBT()))
            plt.xlabel("SheefForce/Tension")
            plt.ylabel("dX/Length")
            plt.errorbar(f,output,yerr=uoutput,fmt='.')
            plt.savefig("dXLengthvsSheefForceTension.png")

        def plotBrownianMotion(self,savePlot = True):

            print("plotting general Motion-----------------")
            if savePlot == True:
                for i in range(len(self.polymers)):
                    for j in range(len(self.polymers[i])):
                        plt.clf()
                        generalProfileWidths = self.polymers[i][j].getWidths()
                        #y = range(len(generalProfileWidths),0,-1)
                        y=[]
                        x=[]
                        mu = self.polymers[i][j].getMean()
                        x = mu[0]
                        y = mu[1]

                        xerr = generalProfileWidths[0]
                        x = [0]*len(x)
                        plt.plot(y,xerr)
                        #plt.errorbar(x,y,xerr=xerr,fmt='o')
                        #plt.set_xlim(-1,1)

                        #t = numpy.linspace(start=x[0]-1, stop=x[-1]+1, num=200)
                        #for k in range(len(x)):

                            #plt.plot(t, self.gaussian(t, x[k] , xerr[k]) + y[k],color='k',linewidth = 1)
                        print(self.name + "_F_" + str(self.forceValues[i]) + "_N_" + str(j) +"_widths.png")

                        plt.savefig(self.name + "_F_" + str(self.forceValues[i]) + "_N_" + str(j) +"_Bwidths.png")
                        #os.system("mv " + self.name + "_p" + str(i)+ ".png " + saveLocation)

        def plotPositionProbabilityData(self,rez=None,Interval=0.5,name = "foo",location=""):
            print("plotting height map")
            frames = []
            if rez == None:
                rez = [100,100]




            indexrange = int(self.parameters.getRunLength() / self.parameters.getProbePeriod())
            self.forceValues = self.getForceRange()

            print(indexrange)
            print(self.forceValues)
            for m in range(len(self.forceValues)):


                p_t = []
                p_t.append([])
                p_t.append([])
                x = numpy.array([])
                y = numpy.array([])
                # for particle in range(len(self.gsd_data[0].particles.position)):
                #     print("F:" + str(round(self.forceValues[m],1)) + ", P:" + str(particle))
                for i in range(int(indexrange*(m + self.interval)),indexrange*(m+1)):

                    x = numpy.append(x,self.gsd_data[i].particles.position[:,0])
                    y = numpy.append(y,self.gsd_data[i].particles.position[:,1])
                    print(i/indexrange)




                #Plot heatmap
                plt.clf()

                plt.title('F:' + str(round(self.forceValues[m],3)) + " probaility map")
                plt.ylabel('y')
                plt.xlabel('x')
                plt.hist2d(x,y,bins=(rez[0],rez[1]),density=True)
                name= "F_" + str(round(self.forceValues[m],3))
                writename = "ProbabilityMap" +name
                plt.savefig(writename + ".png")
                new_frame = Image.open(writename+".png")
                frames.append(new_frame)

            frames[0].save("probabilityMapAnimated.gif", format='GIF',append_images=frames[1:],save_all=True,duration=300,loop=0)
        def plotTilt(self):
            print("plotting tilt")
            tilt = []
            dxdy = []
            sheerForce = []
            for i in range(len(self.forceValues)):
                sheerForce.append(self.forceValues[i]/self.parameters.getPullForce())
                dxdy.append(0)
                for j in range(len(self.polymers[i])):
                    dxdy[-1] += self.polymers[i][j].getdxdyoutput()
                dxdy[-1] = dxdy[-1]/len(self.polymers[i])
            print(sheerForce)
            plt.clf()
            a = numpy.array(sheerForce)
            d = numpy.array(dxdy)

            m,b = numpy.polyfit(a,d,1)


            plt.title('dXdY vs SheerForce/Pull: Ft:' + str(self.parameters.getPullForce()) + ' kbT=' + str(self.parameters.getKBT()))
            plt.xlabel("sheerforce/T")
            plt.ylabel("dx/dy")
            plt.plot(a,m*a+b)
            plt.plot(a,d,'.')
            plt.legend(["slope:"+str(m)])
            plt.show()
            plt.savefig("forceVstilt.png")





class PolymerObject():
    def __init__(self,data,L = 0):
        self.L = L
        self.particle = self.sortData(data)
        self.numParticles = len(self.particle[0][0])
        self.numTime = len(self.particle[0])


        self.generateWidths()


    def getParticleData(self,id):
        return [self.particle[0][id],self.particle[1][id]]

    def getParticleWidth(self,id):
        return self.particleWidth[id]

    def generateWidths(self):
        self.particleWidth = []
        self.particleMean = []
        self.particleStdMean = []

        self.particleWidth = [ self.stdDevX(self.particle[0])  , self.stdDevY(self.particle[1]) ]


        self.particleMean  = [ self.mean(self.particle[0]) , self.mean(self.particle[1]) ]

        #self.particleStdMean.append( [self.particleWidth[-1][0]/math.sqrt(len(self.particle[0][i])),self.particleWidth[-1][1]/len(self.particle[1][i])])

        # for i in range(1,len(self.particleMean)):
        #     # check for dx with particles closer to the center # TODO
        #     if ((self.particleMean[i][0] - self.particleMean[0][0]) < -self.L/2):
        #         self.particleMean[i][0] = self.particleMean[i][0] + self.L
        #     if ((self.particleMean[i][0] - self.particleMean[0][0]) > self.L/2):
        #         self.particleMean[i][0] = self.particleMean[i][0] - self.L
        a = int( 1 / 3 * len(self.particleMean[0]) )
        b = int( 2 / 3 * len(self.particleMean[0]) )
        self.dx = (self.particleMean[0][b] - self.particleMean[0][a])
        self.dy = (self.particleMean[1][b] - self.particleMean[1][a])
        print(type(self.particleMean[0][b]))

        self.calulateLength(a,b)


    def getParticleStdMean(self):
        return self.particleStdMean
    def getDx(self):
        return self.dx
    def getDy(self):
        return self.dy
    def getdxdyoutput(self):
        return self.dx/self.dy
    def mean(self, arr):
        mn = []
        for i in range(len(arr[0])):
            t = 0
            for j in range(len(arr)):
                t += arr[j][i]
            mn.append(t/len(arr))

        return mn

    def calulateLength(self,a,b):
        self.polymerLength = 0


        print(len(self.particleMean));
        x = self.particleMean[0][b] - self.particleMean[0][a]
        y = self.particleMean[1][b] - self.particleMean[1][a]

        self.polymerLength = math.sqrt(x**2 + y**2)
        return self.polymerLength
    def getLength(self):
        return self.polymerLength
    def sortData(self,arr):


        for i in range(len(arr[0][0])):
            for j in range(1,len(arr[0])):
                if arr[0][j][i] - arr[0][j][i-1] < -self.L/2:
                    arr[0][j][i] = (arr[0][j][i] + self.L)
                elif arr[0][j][i] - arr[0][j][i-1] > self.L/2:
                    arr[0][j][i] = (arr[0][j][i] - self.L)
        # for k in range(len(arr)):
        #     for i in range(len(arr[k])):
        #         for j in range(1,len(arr[k][i])):
        #             assert(arr[k][i][j] - arr[k][i][j-1] > -self.L/2)
        #             assert(arr[k][i][j] - arr[k][i][j-1] < self.L/2)

        return arr


    def calculateMeanX(self,arr,particle):
        mean = 0
        for i in range(self.numTime):
            mean += arr[i][particle]
        return mean/len(arr)


    def calculateMeanY(self,arr):
        mean = 0
        for i in range(len(arr)):
            mean += arr[i][particle]
        return mean/len(arr)

    def getMean(self):
        return self.particleMean
    def getWidths(self):
        return self.particleWidth

    def getAngle(self, deg = False):
        if deg:
            return numpy.arccos(self.dx / self.polymerLength)
        else:
            return self.dx / self.polymerLength

    def stdDevX(self,data):
        stdx = []
        for j in range(len(data[0])):
            std = 0
            mean = self.calculateMeanX(data,j)
            difference = 0
            for i in range(len(data)):
                difference += (data[i][j] - mean)**2

            stdx.append(math.sqrt(1 / (len(data) - 1) * difference))

        return stdx
    def stdDevY(self,data):
        stdx = []
        for j in range(len(data[0])):
            std = 0
            mean = self.calculateMeanX(data,j)
            difference = 0
            for i in range(len(data)):
                difference += (data[i][j] - mean)**2

            stdx.append(math.sqrt(1 / (len(data) - 1) * difference))
        return stdx
