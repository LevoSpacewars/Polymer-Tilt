import numpy as np
import matplotlib
matplotlib.use("agg")
from matplotlib import pyplot as plt
import pandas
import matplotlib.cm as cmx
from scipy import stats
import sys
import os
import glob
import subprocess
import matplotlib.backends.backend_pdf
class PolymerSimulationParameters():
    def __init__(self,sheerforcerange=[0,0],df=0,length=0,lines=0,rez=0,K=0,l_0=0,l_max=0,pull=0,amplitude=0,gamma=0,kbT=0,dt=0,probePeriod=0,runLength=0,boxdimx=0,boxdimy=0,integraor = "brownian", direction = "forward"):
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
        return self.direction;


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
        text.write("Direction=" +       str(self.getRunDirection()))
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

        file.close()

COMAVG = []
THETA = []

def makeprofiledist(rawfilepath,dirpath):
    parameter = PolymerSimulationParameters()
    parameter.loadParameters(dirpath + "/_simulation_parameters.txt")
    plength = parameter.getLength()
    npoly = parameter.getNumberChains()
    heatmap = pandas.read_csv(rawfilepath, header=1)
    name = rawfilepath.split(':')[1][0:-4]
    comdata = pandas.read_csv(dirpath + "/COM:" + (name) + ".txt", header=None)
    print(name)
    datauw = np.array(heatmap)
    t_l = int(heatmap.shape[0]/(plength * npoly))
    avgx = []
    avgy = []
    uncx = []

    map = np.array(heatmap)
    x = []
    y = []
    for i in range(len(map)):
        x.append(map[i][0])
        y.append(map[i][1])
    fig = plt.figure()
    bins,a,b,d= plt.hist2d(x,y,bins=(100,100))
    plt.title("Theta: " +str(name))
    pdf.savefig(fig)
    plt.close()

    fig = plt.figure()
    for i in range(plength):
        a = 0
        b = 0
        n = t_l * npoly
        for ii in range(t_l):
            offset = plength * npoly * ii + i
            for iii in range(npoly):
                particle_index = offset + iii*plength
                a += datauw[particle_index][0]
                b += datauw[particle_index][1]
        avgx.append(abs(a/n))
        avgy.append(abs(b/n))

    plt.plot(avgy,avgx)
    plt.yscale('log')
    pdf.savefig(fig)
    plt.close()

    fig = plt.figure()
    plt.plot(avgx,avgy)
    plt.savefig(fig)
    

            






    xcom = []
    ycom = []
    xa = 0
    print(comdata.shape)
    fig = plt.figure()
    ycom.append(0)
    for i in range(1,(comdata.shape[1]-1)):
        xa += pow(comdata[i] - comdata[i-1],2)
        xcom.append(comdata[i])
    COMAVG.append(xa)
    THETA.append(float(name))

    plt.plot(range(len(xcom)),xcom)
    plt.title(" '$(X_{com})$'")
    plt.xlabel(" time ")
    plt.ylabel("x")
    pdf.savefig(fig)
    plt.close()





def makeCOMgraph():
    fig = plt.figure()
    plt.plot(THETA,COMAVG)
    pdf.savefig(fig)
    plt.close()

def makeTILTgraph(directory):



print("RunDirectory:" + sys.argv[1])

runDir = sys.argv[1]

files = os.listdir(runDir)
runDirs = []
for obj in files:
    if(os.path.isdir(obj)):
        runDirs.append(obj)
#print(runDirs)
for i in range(len(runDirs)):
    pdf = matplotlib.backends.backend_pdf.PdfPages(runDirs[i] + "/ProfileDensity.pdf")
    files = glob.glob(runDirs[i] + '/DensityData_raw*.txt')
    files_ordered = []
    order = []
    org_order=[]
    for j in range(len(files)):
        name = float(files[j].split(':')[1][0:-4])
        files_ordered.append(name)
    org_order=files_ordered
    files_ordered = sorted(files_ordered, key = lambda x:float(x))
    for j in range(len(files)):
        order.append(org_order.index(files_ordered[j]))
    
    for j in range(len(order)):
        print(files_ordered[j])
    print(files[j])
    print(files[order[j]])
    print("\n\n")

    for j in range(len(files)):
        makeprofiledist(files[order[j]],runDirs[i])
    makeCOMgraph()
    makeTILTgraph(runDirs[i])
    pdf.close()
   

