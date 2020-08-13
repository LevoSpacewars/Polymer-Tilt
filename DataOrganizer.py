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
    def constructPolymerObjects(self,fileLocation,parameters,interval=0.75):

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
                for k in range(int(indexrange*(m + interval)),indexrange*(m+1)):


                    x.append([])
                    y.append([])

                    for j in range(int(parameters.getLength())):

                        x[-1].append(gsd_data[k].particles.position[i*parameters.getLength() + j,0])

                        y[-1].append(gsd_data[k].particles.position[i*parameters.getLength() + j,1])


                polymers[-1].append(Simulations.PolymerObject([x,y],L=parameters.getNumberChains()))
                print("polymer construction time: " + str(time.perf_counter() - timer))

        print("done Constructing Polymers")
        return polymers
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


                print("polymer construction time: " + str(time.perf_counter() - timer))

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
            self.base = False
            if parameters.getNumberChains() == 1:
                self.base = True


    def isChainLength(self, length: int) -> bool:
        return (self.parameters.getLength() == length)
    def isTemp(self, temp: float) -> bool:
        return (self.parameters.getKBT() == temp)
    def getTemp(self):
        return self.parameters.getKBT()
    def calcForceRange(self):
        forcerange=[]
        for i in range(self.parameters.getDf()):
            forcerange.append(i/self.parameters.getDf() * (self.parameters.getSheerForceRange()[1] - self.parameters.getSheerForceRange()[0]) + self.parameters.getSheerForceRange()[0])
        return forcerange
    def readDataFile(self,fileName):
        file = open(fileName, 'r')
        print(file)
        lines = file.readlines()
        temp = []
        for i in range(len(lines)):
            temp.append([])
            obj = lines[i].rstrip("\n").split(',')
            print(obj[0])
            for i in range(1,len(obj)):
                print(obj[i])
                temp[-1].append(float(obj[i]))
        self.dx         = temp[0]
        self.udx        = temp[1]
        self.length     = temp[2]
        self.ulength    = temp[3]
        self.output     = temp[4]
        self.uoutput    = temp[5]
        self.dy         = temp[6]
    def getChainLength(self):
        return self.parameters.getLength()
    def getDx(self):
        return self.dx,self.udx
    def getLength(self):
        return self.length,self.ulength
    def getOutput(self):
        return self.output, self.uoutput
    def getdxdyoutput(self):
        out = []
        for i in range(len(self.dy)):
            out.append(self.dx[i]/self.dy[i])
        return out

    def getForceRange(self, normalized = False):
        if normalized == True:
            n_forceRange = []
            for i in range(len(self.forceRange)):
                n_forceRange.append(self.forceRange[i]/self.getTension())
            return n_forceRange
        return self.forceRange
    def getTension(self):
        return self.parameters.getPullForce()
    def getPeriodicAmplitude(self):
        return self.parameters.getPeriodicAmplitude();
class GlobalDataPlotter(object):
    # This class should handle EVERYTHING about plotting between data plots

    def __init__(self):
        # should find and store all of the parameter files
        self.dataHandlers = []

        files = glob.glob('**/data.txt',recursive = True)
        for i in range(len(files)):
            dir = files[i].split('/')[0]
            print(dir)
            parameters = Simulations.PolymerSimulationParameters()
            parameters.loadParameters(dir + "/_simulation_parameters.txt")
            self.dataHandlers.append(RunDataHandler(fileName = files[i],parameters = parameters))
        print(len(self.dataHandlers));

    def sortPolymerNumbers(self):
        length = []
        handlers = []
        temp    = []


        length.append(self.dataHandlers[0].getChainLength())
        temp.append(self.dataHandlers[0].getTemp())
        handlers.append([])
        handlers[0].append(self.dataHandlers[0])
        for i in range(1,len(self.dataHandlers)):
            for j in range(len(length)):

                if self.dataHandlers[i].getChainLength() not in length:
                    if self.dataHandlers[i].getTemp() not in temp:
                        temp.append(self.dataHandlers[i].getTemp())
                    length.append(self.dataHandlers[i].getChainLength())
                    handlers.append([])
                    handlers[-1].append(self.dataHandlers[i])
                else:
                    handlers[length.index(self.dataHandlers[i].getChainLength())].append(self.dataHandlers[i])
                    if self.dataHandlers[i].getTemp() not in temp:
                        temp.append(self.dataHandlers[i].getTemp())

        return handlers, length, temp;
    def createColorPallet(self, array):
        num = len(array)
        conv = 1/num
        RGB = []

        for i in range(num):
            x = math.pi/2 * conv * i
            RGB.append([abs(math.sin(x)),abs(math.cos(x)),abs(math.cos(x))])

        print(RGB)
        return RGB
    def createLegendT(self, temp):
        legend = []
        for i in range(len(temp)):
            legend.append("Temperature: " + str(temp[i]))
        return legend

    def createLegendL(self, length):
        legend = []
        for i in range(len(length)):
            legend.append("ChainLength: " + str(length[i]))
        return legend
    def PlotTiltvsForce(self,Tension = 20,kt=1,Dim3=False,base=True):
        from matplotlib.backends.backend_pdf import PdfPages
        from mpl_toolkits.mplot3d import Axes3D
        from matplotlib import pyplot as plt

        sortedDataHandler, length, temp = self.sortPolymerNumbers()
        lengthassignment = []
        colorPallet = self.createColorPallet(length)
        x   = []
        y   = []
        uy  = []

        if base == True :


            bx = []
            by = []
            buy= []

            for j in range(len(self.dataHandlers)):


                if(self.dataHandlers[j].getTension() == Tension) and (self.dataHandlers[j].isTemp(kt) and self.dataHandlers[j].base == False):
                    x.append([])
                    y.append([])
                    uy.append([])

                    bx.append([])
                    by.append([])
                    buy.append([])

                    for p in range(len(self.dataHandlers)):
                        if(self.dataHandlers[p].getTension() == Tension) and self.dataHandlers[p].isTemp(kt) and self.dataHandlers[p].base == True and self.dataHandlers[p].getChainLength() == self.dataHandlers[j].getChainLength() and self.dataHandlers[p].getPeriodicAmplitude() == self.dataHandlers[j].getPeriodicAmplitude():
                            for k in range(len(self.dataHandlers[p].getForceRange(normalized=True))):
                                bx[-1].append(self.dataHandlers[p].getForceRange(normalized=True)[k])
                                by[-1].append(self.dataHandlers[p].getOutput()[0][k])
                                buy[-1].append(self.dataHandlers[p].getOutput()[1][k])



                    for k in range(len(self.dataHandlers[j].getForceRange(normalized=True))):
                        x[-1].append(self.dataHandlers[j].getForceRange(normalized=True)[k])
                        y[-1].append(self.dataHandlers[j].getOutput()[0][k])
                        uy[-1].append(self.dataHandlers[j].getOutput()[1][k])
                    for i in range(len(length)):
                        if(length[i] == self.dataHandlers[j].getChainLength()):

                            lengthassignment.append(i);
                            plt.errorbar(x[-1],y[-1],yerr=uy[-1],fmt='.',c=colorPallet[i])
            plt.legend(self.createLegendL(length))
            plt.title('Polymer dx/dr vs sheerForce/Tension: T = 20, kbT= ' + str(kt));


            fig, axs = plt.subplots(len(x))

            fig.suptitle('Polymer dx/dr vs sheerForce/Tension: T = 20, kbT= ' + str(kt))
            if len(x) <= 1:
                axs.errorbar(x[i], y[i],yerr=uy[i],fmt='.',color = colorPallet[lengthassignment[i]])
                axs.errorbar(bx[i], by[i],yerr=buy[i],fmt='.',color = 'red')

                axs.legend([self.createLegendL(length)[lengthassignment[i]],"BASE Single polymer"])
                plt.text(0.5, 1,'Amplitude=4',
                     horizontalalignment='center',
                     verticalalignment='top',
                     transform = axs.transAxes)
            else:

                for i in range(len(x)):
                    axs[i].errorbar(x[i], y[i],yerr=uy[i],fmt='.',color = colorPallet[lengthassignment[i]])
                    axs[i].errorbar(bx[i], by[i],yerr=buy[i],fmt='.',color = 'red')

                    leg1 = axs[i].legend([self.createLegendL(length)[lengthassignment[i]],"BASE Single polymer"])



            fig.text(0.5, 0.04, 'Sheer/Tension Force', ha='center', va='center')
            fig.text(0.06, 0.5, 'Dx/Length', ha='center', va='center', rotation='vertical')
            plt.savefig('test.pdf')
            plt.show()



        else:
            for j in range(len(self.dataHandlers)):


                if(self.dataHandlers[j].getTension() == Tension) and (self.dataHandlers[j].isTemp(kt)):
                    x.append([])
                    y.append([])
                    uy.append([])
                    print("chain Lenght:")
                    for k in range(len(self.dataHandlers[j].getForceRange(normalized=True))):
                        x[-1].append(self.dataHandlers[j].getForceRange(normalized=True)[k])
                        y[-1].append(self.dataHandlers[j].getOutput()[0][k])
                        uy[-1].append(self.dataHandlers[j].getOutput()[1][k])
                    for i in range(len(length)):
                        if(length[i] == self.dataHandlers[j].getChainLength()):

                            lengthassignment.append(i);
                            plt.errorbar(x[-1],y[-1],yerr=uy[-1],fmt='.',c=colorPallet[i])
            plt.legend(self.createLegendL(length))
            #plt.show()

            fig, axs = plt.subplots(len(x))

            fig.suptitle('Polymer dx/dr vs sheerForce/Tension: T = 20, kbT= ' + str(kt))
            for i in range(len(x)):
                axs[i].errorbar(x[i], y[i],yerr=uy[i],fmt='.',color = colorPallet[lengthassignment[i]])

                axs[i].legend([self.createLegendL(length)[lengthassignment[i]]])


            fig.text(0.5, 0.04, 'Sheer/Tension Force', ha='center', va='center')
            fig.text(0.06, 0.5, 'Dx/Length', ha='center', va='center', rotation='vertical')
            plt.savefig('test.pdf')
            plt.show()

    def PlotTiltvsForce_CL(self,Tension = 10,CL=200,Dim3=False,kt=1):
        from matplotlib.backends.backend_pdf import PdfPages
        from mpl_toolkits.mplot3d import Axes3D
        from matplotlib import pyplot as plt


        print("---------------")
        sortedDataHandler, length, temp = self.sortPolymerNumbers()
        colorPallet = self.createColorPallet(temp)
        print(temp);

        colorassignment = []
        x   = []
        y   = []
        uy  = []

        for j in range(len(self.dataHandlers)):

            if(self.dataHandlers[j].getTension() == Tension) and (self.dataHandlers[j].isChainLength(CL)):
                x.append([])
                y.append([])
                uy.append([])
                print("chain Lenght:")
                for k in range(len(self.dataHandlers[j].getForceRange(normalized=True))):
                    x[-1].append(self.dataHandlers[j].getForceRange(normalized=True)[k])
                    y[-1].append(self.dataHandlers[j].getOutput()[0][k])
                    uy[-1].append(self.dataHandlers[j].getOutput()[1][k])
                for i in range(len(temp)):
                    if(temp[i] == self.dataHandlers[j].getTemp()):

                        colorassignment.append(i);
                        plt.errorbar(x[-1],y[-1],yerr=uy[-1],fmt='.',c=colorPallet[i])
                        print(colorPallet[i])
        plt.legend(self.createLegendT(temp))
        plt.title('Polymer dx/dr vs sheerForce/Tension: T = 20, kbT= ' + str(kt));

        fig, axs = plt.subplots(len(x))

        fig.suptitle('Polymer dx/dr vs sheerForce/Tension: T = 20, kbT= ' + str(kt))
        if(len(x) > 1):
            for i in range(len(x)):
                axs[i].errorbar(x[i], y[i],yerr=uy[i],fmt='.',color = colorPallet[colorassignment[i]])

                axs[i].legend([self.createLegendT(temp)[colorassignment[i]]])
        else:
            axs.errorbar(x[0], y[0],yerr=uy[0],fmt='.',color = colorPallet[colorassignment[0]])

            axs.legend([self.createLegendT(temp)[colorassignment[i]]])


        fig.text(0.5, 0.04, 'Sheer/Tension Force', ha='center', va='center')
        fig.text(0.06, 0.5, 'Dx/Length', ha='center', va='center', rotation='vertical')
        plt.savefig('test.pdf')
        plt.show()
    def compareBases_CL(self,Tension = 2,CL=200,Dim3=False,kt=0,Amplitude = 0):
        from matplotlib import pyplot as plt
        import numpy
        bx = []
        by = []
        buy = []


        print(len(self.dataHandlers))

        for p in range(len(self.dataHandlers)):
            if(self.dataHandlers[p].getTension() == Tension) and self.dataHandlers[p].isTemp(kt) and self.dataHandlers[p].base == True and self.dataHandlers[p].getPeriodicAmplitude() == Amplitude:
                bx.append([])
                by.append([])
                buy.append([])
                for k in range(len(self.dataHandlers[p].getForceRange(normalized=True))):
                    bx[-1].append(self.dataHandlers[p].getForceRange(normalized=True)[k])
                    by[-1].append(self.dataHandlers[p].getdxdyoutput()[k])
        m=[]
        for i in range(len(bx)):
            print('here')

            plt.plot(bx[i],by[i],'o')
            x = numpy.array(bx[i])
            y = numpy.array(by[i])
            ma,b = numpy.polyfit(x,y,1)
            m.append(ma)
            plt.plot(x,ma*x + b)
        plt.ylabel("Dx/Dy")
        plt.xlabel("Fs/T")
        plt.title("L:100,Ft:2,kT:0,SinglePolymer,A:0")
        plt.legend(["measured","slope:" + str(m[0])])
        plt.show();



# compile = GlobalDataManager()
# compile.compileNew()

plotter = GlobalDataPlotter()
plotter.PlotTiltvsForce_CL(Tension=20)
