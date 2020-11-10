import glob
import Simulations
import time
import os
import gsd.hoomd
import math
from matplotlib.backends.backend_pdf import PdfPages
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
            dir = files[i].split('/')[0]
            print(dir)
            parameters = Simulations.PolymerSimulationParameters()
            parameters.loadParameters(dir + "/_simulation_parameters.txt")
            self.dataHandlers.append(RunDataHandler(fileName = files[i],parameters = parameters))
        print(len(self.dataHandlers));

    def sortPolymerNumbers(self):
        temp    = []
        amplitude = []
        length = []

        handlers = []
        thandlers = []
        ahandlers = []
        


        length.append(self.dataHandlers[0].getChainLength())
        temp.append(self.dataHandlers[0].getTemp())
        amplitude.append(self.dataHandlers[0].getPeriodicAmplitude())
        
        handlers.append([])
        thandlers.append([])
        ahandlers.append([])
        
        handlers[0].append(self.dataHandlers[0])
        thandlers[0].append(self.dataHandlers[0])
        ahandlers[0].append(self.dataHandlers[0])
        
        for i in range(1,len(self.dataHandlers)):
            for j in range(len(length)):



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



                if self.dataHandlers[i].getChainLength() not in length:                
                    length.append(self.dataHandlers[i].getChainLength())
                    handlers.append([])
                    handlers[-1].append(self.dataHandlers[i])
                else:
                    handlers[length.index(self.dataHandlers[i].getChainLength())].append(self.dataHandlers[i])
                    



        return handlers, thandlers, ahandlers, length, temp, amplitude;
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
    def compareBases_CL(self,Tension = 20,CL=200,Dim3=False,kt=1,Amplitude = 3):
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
                    by[-1].append(self.dataHandlers[p].getOutput()[0][k])
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
        plt.title("L:" + str(CL) + ",Ft:" + str(Tension) + ",kT:"+str(kt)+",SinglePolymer,A:"+str(Amplitude))
        plt.legend(["measured","slope:" + str(m[0])])
        plt.show();
    def PlotTiltvsForceX_CL(self,Tension = 10,CL=200,Dim3=False,kt=1):
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
                    y[-1].append(self.dataHandlers[j].getDx()[0][k])
                    uy[-1].append(self.dataHandlers[j].getDx()[1][k])
                for i in range(len(temp)):
                    if(temp[i] == self.dataHandlers[j].getTemp()):

                        colorassignment.append(i);
                        plt.errorbar(x[-1],y[-1],yerr=uy[-1],fmt='.',c=colorPallet[i])
                        print(colorPallet[i])
        plt.legend(self.createLegendT(temp))
        plt.title('Polymer dx/r vs sheerForce/Tension: T = 20, kbT= ' + str(kt));

        fig, axs = plt.subplots(len(x))

        fig.suptitle('Polymer dx/r vs sheerForce/Tension: T = 20, kbT= ' + str(kt))
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


    def plotdtiltvsdforce_CL(self,Tension = 10,CL=200,Dim3=False,kt=1):

        from matplotlib.backends.backend_pdf import PdfPages
        from mpl_toolkits.mplot3d import Axes3D
        from matplotlib import pyplot as plt


        print("---------------")
        sortedDataHandler, length, temp = self.sortPolymerNumbers()
        colorPallet = self.createColorPallet(temp)
        print(temp);

        colorassignment = []
        df = []
        do = []
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
                    y[-1].append(self.dataHandlers[j].getDx()[0][k])
                    uy[-1].append(self.dataHandlers[j].getDx()[1][k])
                for i in range(len(temp)):
                    if(temp[i] == self.dataHandlers[j].getTemp()):
                        colorassignment.append(i);
                        #plt.errorbar(x[-1],y[-1],yerr=uy[-1],fmt='.',c=colorPallet[i])
                        print(colorPallet[i])

                df.append([])
                do.append([])

                for k in range(len(x[-1])-1):
                    df[-1].append(x[-1][k+1] - x[-1][k])
                    do[-1].append(y[-1][k+1] - y[-1][k])
                    do[-1][-1] = do[-1][-1]/df[-1][-1]

        plt.legend(self.createLegendT(temp))
        plt.title('Polymer dx/r vs sheerForce/Tension: T = 20, kbT= ' + str(kt))

        fig, axs = plt.subplots(len(x))

        if(len(x) > 1):
            for i in range(len(x)):
                axs[i].errorbar(x[i][:-1], do[i],fmt='.',color = colorPallet[colorassignment[i]])

                axs[i].legend([self.createLegendT(temp)[colorassignment[i]]])
        else:
            axs.errorbar(df[0], do[0],yerr=uy[0],fmt='.',color = colorPallet[colorassignment[0]])

            axs.legend([self.createLegendT(temp)[colorassignment[i]]])
        plt.show()

    

    def GraphEverythingToOrganizedPDF(self):
        from matplotlib import pyplot as plt
        slength, stemp, samplitude, list_l, list_t, list_a = self.sortPolymerNumbers()
        import matplotlib.patches as mpatches
        #By Tempterature
        # first by identical parameters
        tfigures = []
        ttitles = []
        list_colors = ["aqua","black","chocolate", "blue", "green","grey","red"]
        with PdfPages("exportTest.pdf") as export_pdf:
            for i in range(len(list_a)):
                for j in range(len(list_l)):
                    title = "Varied Temperature with Length:" + str(list_l[j]) + ", Amplitude:" + str(list_a[i])
                    
                    ttitles.append(title)
                    fig = (plt.figure(title))
                    tfigures.append(fig)
                    plt.title(title)
                    plt.xlabel("sheerforce/tension")
                    plt.ylabel("dx/length")
                    lc=[]
                    legendlabel=[]
                    for k in range(len(stemp)):
                        
                        for l in range(len(stemp[k])):
                            if list_a[i] == stemp[k][l].getPeriodicAmplitude() and list_l[j] == stemp[k][l].getChainLength():
                                x = stemp[k][l].getForceRange()
                                y,uy = stemp[k][l].getOutput()
                                color =""
                                if list_colors[list_t.index(stemp[k][l].getTemp())] not in lc:
                                    lc.append(list_colors[list_t.index(stemp[k][l].getTemp())])
                                    legendlabel.append("MP: KbT="+str(stemp[k][l].getTemp()))
                                color = list_colors[list_t.index(stemp[k][l].getTemp())]
                                if stemp[k][l].base is True:
                                    if "red" not in lc:
                                        lc.append("red")
                                        legendlabel.append("SP: KbT="+str(stemp[k][l].getTemp()))
                                        color = "red"
                                
                                
                                
                                plt.errorbar(x,y,yerr=uy,color=color)
                    patches = []
                    for z in range(len(lc)):
                        patches.append(mpatches.Patch(color=lc[z], label=legendlabel[z]))
                    plt.legend(handles=patches)
                    export_pdf.savefig()
                    plt.close()
                    plt.clf()
            
                    



        
            tfigures.append(plt.figure("VariedTemperature_constant_Merged"))
        #By Length
            lfigures = []
            ltitles = []
            for i in range(len(list_a)):
                    for j in range(len(list_t)):
                        title = "Varied length with Temperature:" + str(list_t[j]) + ", Amplitude:" + str(list_a[i])
                        
                        ttitles.append(title)
                        fig = (plt.figure(title))
                        tfigures.append(fig)
                        plt.title(title)
                        plt.xlabel("sheerforce/tension")
                        plt.ylabel("dx/length")
                        for k in range(len(stemp)):
                            for l in range(len(stemp[k])):
                                if list_a[i] == stemp[k][l].getPeriodicAmplitude() and list_t[j] == stemp[k][l].getTemp():
                                    x = stemp[k][l].getForceRange()
                                    y,uy = stemp[k][l].getOutput()
                                    color = list_colors[list_l.index(stemp[k][l].getChainLength())]
                                    
                                    plt.errorbar(x,y,yerr=uy,color=color)
                        export_pdf.savefig(figure=fig)
                        plt.close()
                        plt.clf()
        # lfigures.append(plt.figure("VariedLength_constant_Temperature:" + str()))
        # lfigures.append(plt.figure("VariedLength_constant_Merged"))
        #By Amplitude 
            lfigures = []
            ltitles = []
            for i in range(len(list_t)):
                    for j in range(len(list_l)):
                        title = "Varied Amplitude with Length:" + str(list_l[j]) + ", Temperature:" + str(list_t[i])
                        
                        ttitles.append(title)
                        fig = (plt.figure(title))
                        tfigures.append(fig)
                        plt.title(title)
                        plt.xlabel("sheerforce/tension")
                        plt.ylabel("dx/length")
                        for k in range(len(stemp)):
                            for l in range(len(stemp[k])):
                                if list_t[i] == stemp[k][l].getTemp() and list_l[j] == stemp[k][l].getChainLength():
                                    x = stemp[k][l].getForceRange()
                                    y,uy = stemp[k][l].getOutput()
                                    color = list_colors[list_t.index(stemp[k][l].getTemp())]
                                    
                                    plt.errorbar(x,y,yerr=uy,color=color)
                        export_pdf.savefig(figure=fig)
                        plt.close()
                        plt.clf()

        
        #All Together

        #Export to pdf





# compile = GlobalDataManager()
# compile.compileNew()

plotter = GlobalDataPlotter()
# plotter.PlotTiltvsForce(Tension=20)
# plotter.compareBases_CL(Tension = 20,Amplitude = 0.3/0.1 * 1,kt=1)
plotter.GraphEverythingToOrganizedPDF()
