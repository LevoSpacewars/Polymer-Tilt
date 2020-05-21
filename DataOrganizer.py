import glob
import Simulations
import time
import os
import gsd.hoomd
import math
class GlobalDataManager(object):
    # compile new data into files
    # recompile all data into files
    #this is a helper class for GlobalDataPlotter;
    def __init__(self):
        pass
    def compileNew(self):
        files = glob.glob('**/trajectory.gsd',recursive = True)
        for i in range(len(files)):
            dir = files[i].split('/')[0]

            if (os.path.isfile(dir + "/data.txt")) == False:
                print(files[i] + ": has no data.txt")
                self.writeData(files[i],dir)


    def recompileAll(self):
        files = glob.glob('**/trajectory.gsd',recursive = True)
        for i in range(len(files)):
            dir = files[i].split('/')[0]

            if (os.path.isfile(dir + "/data.txt")) == False:
                print(files[i] + ": has no data.txt")
                self.writeData(files[i],dir)

            else:
                os.system("rm "+ dir + "/data.txt")
                print(files[i] + ": removed data.txt")
                self.writeData(files[i],dir)
    def recompile(self, filePath):
        dir = filePath.split('/')[0]
        if (os.path.isfile(dir + "/data.txt")) == False:
            print(filePath + ": has no data.txt")
            self.writeData(filePath,dir)

        else:
            os.system("rm "+ dir + "/data.txt")
            print(filePath + ": removed data.txt")
            self.writeData(filePath,dir)

    def writeData(self,path,parentDirectory):
        parameters = self.getSimulationParameters(parentDirectory)
        polymers   = self.constructPolymerObjects(path,parameters)

        file = open(parentDirectory + "/data.txt","w+")

        dx,udx,length,ulength,output,uoutput = self.get_dx_length_output(polymers)
        s_dx = "dx," + self.formatArray(dx) + "\n"
        s_udx = "udx," + self.formatArray(udx) + "\n"
        s_length= "length," + self.formatArray(length) + "\n"
        s_ulength= "ulength," + self.formatArray(ulength) + "\n"
        s_output= "output," + self.formatArray(output) + "\n"
        s_uoutput= "uoutput," + self.formatArray(uoutput) + "\n"
        file.write(s_dx)
        file.write(s_udx)
        file.write(s_length)
        file.write(s_ulength)
        file.write(s_output)
        file.write(s_uoutput)
        file.close()
        #

    def getSimulationParameters(self,dir):
        parameters = Simulations.PolymerSimulationParameters()
        parameters.loadParameters(fileLocation=(dir +"/_simulation_parameters.txt"))
        return parameters
    def constructPolymerObjects(self,fileLocation,parameters,interval=0.0):

        print("constructing polymer Objects")

        forceValues = self.getForceRange(parameters)
        gsd_data = gsd.hoomd.open(fileLocation,'rb')
        timer = 0
        indexrange = int(parameters.getRunLength() / parameters.getProbePeriod())
        polymers = []

        for m in range(len(forceValues)):
            polymers.append([])

            for i in range(int(parameters.getNumberChains())):
                timer = time.perf_counter()

                print(str(m) + "." + str(i))
                x = []
                y = []

                for j in range(int(parameters.getLength())):

                    x.append([])
                    y.append([])
                    for k in range(int(indexrange*(m + interval)),indexrange*(m+1)):

                        x[-1].append(gsd_data[k].particles.position[i*parameters.getLength() + j,0])

                        y[-1].append(gsd_data[k].particles.position[i*parameters.getLength() + j,1])


                polymers[-1].append(Simulations.PolymerObject([x,y],L=parameters.getNumberChains()))
                print("polymer construction time: " + str(time.perf_counter() - timer))

        print("done Constructing Polymers")
        return polymers

    def getForceRange(self, parameters):
        forcerange=[]
        for i in range(parameters.getDf()):
            forcerange.append(i/parameters.getDf() * (parameters.getSheerForceRange()[1] - parameters.getSheerForceRange()[0]) + parameters.getSheerForceRange()[0])
        return forcerange

    def get_dx_length_output(self,polymers):
        length=[]
        ulength=[]

        dx=[]
        udx=[]
        output = []
        uoutput = []
        for m in range(len(polymers)):
            dx.append(0)
            length.append(0)
            udx.append(0)
            ulength.append(0)
            output.append(0)
            uoutput.append(0)
            for j in range(len(polymers[m])):
                dx[-1] += polymers[m][j].getDx()
                length[-1] += polymers[m][j].getLength()
            dx[-1] = dx[-1]/len(polymers[m])
            length[-1] = length[-1]/len(polymers[m])
            output[-1] = dx[-1]/length[-1]
            for j in range(len(polymers[m])):
                udx[-1] += (polymers[m][j].getDx() - dx[-1])**2
                ulength[-1] += (polymers[m][j].getLength() - length[-1])**2
            udx[-1] = math.sqrt( udx[-1] * 1/(len(polymers[m])-1) )
            ulength[-1] = math.sqrt( ulength[-1] * 1/(len(polymers[m])-1) )
            square = ( udx[-1] / dx[-1] )**2 + ( ulength[-1] / length[-1] )**2
            uoutput[-1] = output[-1] * math.sqrt(square)

        return dx,udx,length,ulength,output,uoutput

    def formatArray(self,array):
        fmt =""
        for i in range(len(array)-1):
            fmt += str(array[i]) + ","
        fmt += str(array[-1])
        return fmt
class RunDataHandler(object):
    def __init__(self, fileName = "",allData=[],parameters = None):
        self.parameters = parameters
        if fileName is "":
            self.dx         = allData[0]
            self.udx        = allData[1]
            self.length     = allData[2]
            self.ulength    = allData[3]
            self.output     = allData[4]
            self.uoutput    = allData[5]
            self.forceRange = self.calcForceRange()
        else:
            self.readDataFile(fileName)
            self.forceRange = self.calcForceRange()

    def isChainLength(self, length: int) -> bool:
        return (self.parameters.getLength() == length)
    def isTemp(self, temp: float) -> bool:
        return (self.parameters.getKBT())
    def calcForceRange(self):
        forcerange=[]
        for i in range(self.parameters.getDf()):
            forcerange.append(i/self.parameters.getDf() * (self.parameters.getSheerForceRange()[1] - self.parameters.getSheerForceRange()[0]) + self.parameters.getSheerForceRange()[0])
        return forcerange
    def readDataFile(self,fileName):
        file = open(fileName, 'r')
        lines = file.readlines()
        temp = []
        for i in range(len(lines)):
            temp.append([])
            obj = lines[i].rstrip("\n").split(',')
            for i in range(1,len(obj)):
                print(obj[i])
                temp[-1].append(float(obj[i]))
        self.dx         = temp[0]
        self.udx        = temp[1]
        self.length     = temp[2]
        self.ulength    = temp[3]
        self.output     = temp[4]
        self.uoutput    = temp[5]
    def getChainLength(self):
        return self.parameters.getLength()
    def getDx(self):
        return self.dx,self.udx
    def getLength(self):
        return self.length,self.ulength
    def getOutput(self):
        return self.output, self.uoutput
    def getForceRange(self, normalized = False):
        if normalized == True:
            n_forceRange = []
            for i in range(len(self.forceRange)):
                n_forceRange.append(self.forceRange[i]/self.getTension())
            return n_forceRange
        return self.forceRange
    def getTension(self):
        return self.parameters.getPullForce()
class GlobalDataPlotter(object):
    # This class should handle EVERYTHING about plotting between data plots

    def __init__(self):
        # should find and store all of the parameter files
        self.dataHandlers = []

        files = glob.glob('**/data.txt',recursive = True)
        for i in range(len(files)):
            dir = files[i].split('/')[0]
            parameters = Simulations.PolymerSimulationParameters()
            parameters.loadParameters(dir + "/_simulation_parameters.txt")
            self.dataHandlers.append(RunDataHandler(fileName = files[i],parameters = parameters))
    def sortPolymerNumbers(self):
        length = []
        handlers = []


        length.append(self.dataHandlers[0].getChainLength())
        handlers.append([])
        handlers[0].append(self.dataHandlers[0])
        for i in range(1,len(self.dataHandlers)):
            for j in range(len(length)):
                if self.dataHandlers[i].getChainLength() not in length:
                    length.append(self.dataHandlers[i].getChainLength())
                    handlers.append([])
                    handlers[-1].append(self.dataHandlers[i])
                else:
                    handlers[length.index(self.dataHandlers[i].getChainLength())].append(self.dataHandlers[i])
        for i in range(len(handlers)):
            print(len(handlers[i]))
        return handlers, length
    def createColorPallet(self, array):
        num = len(array)
        conv = 1/num
        RGB = []

        for i in range(num):
            x = math.pi * conv * i
            RGB.append([math.sin(x),math.cos(x),math.cos(x)*math.sin(x)])

        print(RGB)
        return RGB
    def createLegend(self, length):
        legend = []
        for i in range(len(length)):
            legend.append("ChainLength: " + str(length[i]))
        return legend
    def PlotTiltvsForce(self,Tension = 20,kt=1,Dim3=False):
        from matplotlib.backends.backend_pdf import PdfPages
        from mpl_toolkits.mplot3d import Axes3D
        from matplotlib import pyplot as plt

        if Dim3 == True:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            for i in range(len(self.dataHandlers)):
                if(self.dataHandlers[i].getTension() == Tension):
                    xs = [self.dataHandlers[i].getChainLength()] * len(self.dataHandlers[i].getForceRange())
                    ys = self.dataHandlers[i].getForceRange(normalized=True)
                    zs = self.dataHandlers[i].getOutput()[0]

                    ax.scatter(xs,ys,zs,'o')
            plt.show()
            pass
        else:
            sortedDataHandler, length = self.sortPolymerNumbers()
            colorPallet = self.createColorPallet(length)

            x=[]
            y=[]
            for i in range(len(length)):
                x.append([])
                y.append([])
                for j in range(len(sortedDataHandler[i])):
                    if(sortedDataHandler[i][j].getTension() == Tension):
                        for k in range(len(sortedDataHandler[i][j].getForceRange(normalized=True))):
                            x[-1].append(sortedDataHandler[i][j].getForceRange(normalized=True)[k])
                            y[-1].append(sortedDataHandler[i][j].getOutput()[0][k])

                plt.plot(x[-1],y[-1],'.',c=colorPallet[i])
            plt.legend(self.createLegend(length))
            #plt.show()

            fig, axs = plt.subplots(len(x))

            fig.suptitle('Vertically stacked ')
            for i in range(len(x)):
                axs[i].plot(x[i], y[i],'.')
                axs[i].legend([self.createLegend(length)[i]])


            fig.text(0.5, 0.04, 'Sheer/Tension Force', ha='center', va='center')
            fig.text(0.06, 0.5, 'Dx/Length', ha='center', va='center', rotation='vertical')
            plt.show()
            plt.savefig('test.pdf')


compile = GlobalDataManager()
compile.compileNew()

plotter = GlobalDataPlotter()
plotter.PlotTiltvsForce()
