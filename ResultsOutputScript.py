import glob
import Simulations
import time
import os
import gsd.hoomd
import math
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.backends.backend_pdf
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches

def X (a,kbt,A,T,theta):
    theta = math.atan(theta)
    print((A/0.2)*math.tanh(1/(a*10*2*math.pi)))
    print(A)


    F = math.sqrt( pow(T,2.0) + pow(T * math.tan(theta),2.0))
    B = 1/(a*10*2*math.pi)
    am = A * 1/(a*10*2*math.pi)

    V = abs(am)/0.2


    output = F*pow(a,2) * pow(1/kbt,2) * V / pow(math.pi,2)
    return output

def Y (a,kbt,T,theta,utheta):
    theta = math.atan(theta)
    utheta = math.atan(utheta)
    F = math.sqrt( T**2 + pow(T * math.tan(theta),2))
    output =F*a*theta * (1/kbt) * 1/math.pi
    uoutput = F*a*utheta * (1/kbt) * 1/math.pi
    return (output, uoutput)


class RunDataHandler(object):
    def __init__(self, fileName = "",allData=None,parameters = None):
        self.parameters = parameters
        if fileName is "":
            self.dx         = allData[0]
            self.udx        = allData[1]
            self.length     = allData[2]
            self.ulength    = allData[3]
            self.output     = allData[4]
            self.uoutput    = allData[5]
            self.dx2        = allData[6]
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
        self.fileName = fileName
        file = open(fileName, 'r')
        print(file)
        lines = file.readlines()
        print(fileName)
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
        self.dx2        = temp[6]
    def getChainLength(self):
        return self.parameters.getLength()
    def getDx(self):
        return self.dx,self.udx
    def getDx2(self):
        return self.dx2
    def getLength(self):
        return self.length,self.ulength
    def getOutput(self):
        return self.output, self.uoutput
    def getdxdyoutput(self):
        out = []
        for i in range(len(self.dy)):
            out.append(self.dx[i]/self.dy[i])
        return out

    def getForceRange(self, normalized = True):
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
            dir = '/'.join(files[i].split('/')[:-1])
            print(dir)
            parameters = Simulations.PolymerSimulationParameters()
            parameters.loadParameters(dir + "/_simulation_parameters.txt")
            self.dataHandlers.append(RunDataHandler(fileName = files[i],parameters = parameters))
        print(len(self.dataHandlers));

    def sortPolymerNumbers(self):
        temp      = []
        amplitude = []
        length    = []
        tension   = []

        handlers  = []
        thandlers = []
        ahandlers = []
        tnhandlers= []


        length.append(self.dataHandlers[0].getChainLength())
        temp.append(self.dataHandlers[0].getTemp())
        amplitude.append(self.dataHandlers[0].getPeriodicAmplitude())
        tension.append(self.dataHandlers[0].getTension())

        handlers.append([])
        thandlers.append([])
        ahandlers.append([])
        tnhandlers.append([])

        handlers[0].append(self.dataHandlers[0])
        thandlers[0].append(self.dataHandlers[0])
        ahandlers[0].append(self.dataHandlers[0])
        tnhandlers[0].append(self.dataHandlers[0])

        for i in range(1,len(self.dataHandlers)):
            #for j in range(len(length)): ?
            if self.dataHandlers[i].getTemp() not in temp:
                    temp.append(self.dataHandlers[i].getTemp())
                    thandlers.append([])
                    thandlers[-1].append(self.dataHandlers[i])
            else:
                thandlers[temp.index(self.dataHandlers[i].getTemp())].append(self.dataHandlers[i])



            if self.dataHandlers[i].getPeriodicAmplitude() not in amplitude:
                    amplitude.append(self.dataHandlers[i].getPeriodicAmplitude())
                    ahandlers.append([])
                    ahandlers[-1].append(self.dataHandlers[i])
            else:
                ahandlers[amplitude.index(self.dataHandlers[i].getPeriodicAmplitude())].append(self.dataHandlers[i])


            if self.dataHandlers[i].getTension() not in tension:
                tension.append(self.dataHandlers[i].getTension())
                tnhandlers.append([])
                tnhandlers[-1].append(self.dataHandlers[i])
            else:
                tnhandlers[tension.index(self.dataHandlers[i].getTension())].append(self.dataHandlers[i])


            if self.dataHandlers[i].getChainLength() not in length:
                length.append(self.dataHandlers[i].getChainLength())
                handlers.append([])
                handlers[-1].append(self.dataHandlers[i])
            else:
                handlers[length.index(self.dataHandlers[i].getChainLength())].append(self.dataHandlers[i])


        return handlers, thandlers, ahandlers, tnhandlers, length, temp, amplitude, tension

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
    def createLegendTn(self, tension):
        legend = []
        for i in range(len(temp)):
            legend.append("Tension: " + str(tension[i]))
        return legend
    def createLegendL(self, length):
        legend = []
        for i in range(len(length)):
            legend.append("ChainLength: " + str(length[i]))
        return legend


    def GraphEverythingImporant(self):
        from matplotlib import pyplot as plt
        slength, stemp, samplitude,stension,  list_l, list_t, list_a, list_tn = self.sortPolymerNumbers()
        import matplotlib.patches as mpatches
        pdf = matplotlib.backends.backend_pdf.PdfPages("exportTestB.pdf")
        list_colors = ["aqua","black","chocolate", "blue", "green","grey","red"]





        for tn in (list_tn):
            for tmp in (list_t):
                for amp in (list_a):
                    render = False
                    fig = plt.figure()
                    lengths = []
                    for i in range(len(slength)):
                        for polymer in slength[i]:
                            if(polymer.getTemp() == tmp and polymer.getTension() == tn and polymer.getPeriodicAmplitude() == amp):

                                color = list_l.index(polymer.getChainLength())
                                x = polymer.getForceRange()
                                y, uy = polymer.getOutput()
                                print(len(x),len(y))
                                plt.errorbar(x,y,yerr=uy,color=list_colors[color])
                                if polymer.getChainLength in lengths:
                                    lengths.append(polymer.getChainLength())
                                render = True



                    title = "dlength, Temperature:" + str(tmp) + ", Amplitude:" + str(abs(round(amp,3))) + ", $F_{y}$:" + str(tn)
                    plt.title(title)
                    plt.ylabel("dx/length")
                    plt.xlabel("sheerforce/tension")
                    labels = self.createLegendL(list_l)
                    patches = []
                    for i in range(len(labels)):
                        patches.append(mpatches.Patch(color=list_colors[i], label=labels[i]))
                    plt.legend(handles=patches)
                    if(render):
                        pdf.savefig(fig)
                    plt.close()


        pdf.close()


        #graph of variable lengths, c temp, c tension, c a



    def PlotResults(self): #jump to line 391 for where I plot the data
        # if you want to see what properties the polymer simulation has go to the RunDataHandler class definition to see what you can access
        import numpy as np
        from matplotlib import pyplot as plt
        import matplotlib.patches as mpatches
        pdf = matplotlib.backends.backend_pdf.PdfPages("method_test.pdf")

        cmap = plt.get_cmap('jet')
        colors = cmap(np.linspace(0, 1.0, len(self.dataHandlers)))
        cindex = 0
        fig = plt.figure(figsize=(15,15))
        # s is the points for the boundary curve generated by Abhijeet
        s = "{{0., 0.}, {0.3, 0.149379}, {0.6, 0.295278}, {0.9, 0.435244}, {1.2, 0.568148}, {1.5, 0.693871}, {1.8, 0.812837}, {2.1, 0.925679}, {2.4, 1.03306}, {2.7, 1.13559}, {3., 1.23382}, {3.3, 1.32821}, {3.6, 1.41916}, {3.9, 1.50701}, {4.2, 1.59204}, {4.5, 1.67452}, {4.8, 1.75466}, {5.1, 1.83264}, {5.4, 1.90864}, {5.7, 1.98279}, {6., 2.05523}, {6.3, 2.12607}, {6.6, 2.19541}, {6.9, 2.26335}, {7.2, 2.32996}, {7.5, 2.39533}, {7.8, 2.45952}, {8.1, 2.52259}, {8.4, 2.58461}, {8.7, 2.64562}, {9., 2.70566}, {9.6,2.82306}, {9.9, 2.88049}, {10.2, 2.93711}, {10.5, 2.99297}, {10.8, 3.04809}, {11.1, 3.10251}, {11.4, 3.15624}, {11.7, 3.20932}, {12., 3.26177}, {12.3, 3.31361}, {12.6, 3.36486}, {12.9, 3.41554}}"
        s = s.replace("{","")
        s = s.replace("}","")
        s = s.split(",")
        rowx = []
        rowy = []
        for i in range(0,len(s),2):
            rowx.append(float(s[i]))
            rowy.append(float(s[i+1]))
        plt.plot(rowx,rowy)




        sloperangedict = {} #contains for sloperange[i] slopeinfo[i] = (theta, tilt), (...), (....)

        # this loops calculates the results values for a run
        for i in self.dataHandlers:

            fig, axs= plt.subplots(2,3)
            fig.set_figheight(10)
            fig.set_figwidth(20)
            thetas = i.getForceRange()
            dx2 = i.getDx2()
            tilt = i.getOutput()[0]
            utilt = i.getOutput()[1]
            tx = 0
            ty = 0
            dx = 0
            dy = 0
            tiltc = 0
            dxc = 0
            un = i.udx.copy()
            un2 = i.udx.copy()
            x = i.getForceRange()
            filename = i.fileName
            print(filename)

            for j in range(len(un2)):
                un[j] = abs(un[j])
                un2[j] *= i.getOutput()[0][j] * 100

            sloperange = np.array(range(1,6))/10
            for k in range(1,len(tilt)):

                if tilt[k] / abs(tilt[0]) > 10 and tiltc == 0:
                    ty = (tilt[k])
                    tx = (thetas[k])
                    tiltc = 1
                if dx2[k] / abs(dx2[0]) > 1.1 and dxc == 0:
                    dy = (dx2[k])
                    dx = (thetas[k])
                    dxc = 1

            linetilt = [] #list of slope details (slope, theta, tilt)
            linedx = []

            for slope in sloperange:

                a = 1
                kbt = i.getTemp()
                A = i.getPeriodicAmplitude()
                T = i.getTension()
                utheta = i.getForceRange()[1] - i.getForceRange()[0]

                print(len(thetas),len(tilt))

                linetilt.append(self.getIntersection(thetas, tilt, slope))
                if(slope not in sloperangedict):
                    sloperangedict[slope] = []
                if linetilt[-1] is not None:
                    tiltc = linetilt[-1][1]
                    x = X(a, kbt, A, T, tiltc)
                    y= Y(a, kbt, T, tiltc, utheta)[0]
                    uy = Y(a, kbt, T, tiltc, utheta)[1]
                    sloperangedict[slope].append((abs(x),y,uy))

                linedx.append(self.getIntersection(thetas,dx2,slope))

            done = False
            while not done:
                done = True
                if None  in linetilt:
                    linetilt.remove(None)
                    done *= False

                if None  in linedx:
                    linedx.remove(None)
                    done *= False





        keys = sloperangedict.keys()
        
        ################################### RESULTS GRAPH START#############################################

        cmap = plt.get_cmap('jet')
        colors = cmap(np.linspace(0, 1.0, 5))

        figure = plt.figure(figsize=(10,10))
        plt.plot(rowx,rowy)
        index = 0



        for point in sloperangedict[0.1]:
            plt.errorbar(point[0],point[1],fmt='.' ,yerr=point[2],color=colors[index])
        index +=1

        pdf.savefig(figure)
        plt.close(figure)


        pdf.close()
        #################################### RESULTS GRAPH END  ############################################

    def writeIntersectionData(self):

        wfile = open("update_points",'w')
        info = []

        for i in self.dataHandlers:
            tilt = i.getOutput()[0]
            slope = 0.1
            thetas = i.getForceRange()

            a = 1
            kbt = i.getTemp()
            A = i.getPeriodicAmplitude()
            T = i.getTension()
            L = i.getChainLength()

            utheta = i.getForceRange()[1] - i.getForceRange()[0]



            linetilt = (self.getIntersection(thetas, tilt, slope))
            if linetilt is not None:
                info.append((a,kbt,A,T,linetilt[1],utheta,L))

        for i in info:
            wfile.write(f"{i[0]}, {i[1]}, {i[2]}, {i[3]}, {i[4]}, {i[5]}, {i[6]}\n")

        wfile.close()





    def getIntersection(self, xd,yd, slope):

        thetas = xd
        tilt = yd
        l2 = [(0,0),(0.2,slope*0.2)]
        for p in range(3,len(thetas)):

            p2 = tilt[p]
            if p2 > slope * thetas[p]:
                return (slope, thetas[p], tilt[p])



        return None

    def genLine(self, slope, end):
        return [[0,end],[0,slope*end]]





    def line_intersection(self,line1, line2):
        xdiff = (line1[0][0] - line1[1][0], line2[0][0] - line2[1][0])
        ydiff = (line1[0][1] - line1[1][1], line2[0][1] - line2[1][1])

        def det(a, b):
            return a[0] * b[1] - a[1] * b[0]
        tn = True
        div = det(xdiff, ydiff)
        if div == 0:
            tn = False

        d = (det(*line1), det(*line2))
        x = det(d, xdiff) / div
        y = det(d, ydiff) / div
        return (tn, x, y)



plotter = GlobalDataPlotter()

plotter.PlotResults()
