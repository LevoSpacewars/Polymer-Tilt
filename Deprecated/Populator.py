import hoomd
import hoomd.md
import Util
import math
from Util import add
class Polymer:
    #class attributes
    classID = Util.ID
    positions   = [];              # holds all the global positions of every particle
    ids         = [];
    bonds       = [];
    def __init__(self, iN, name="polymer"):
        self.classID.name = name
        self.classID.IN = iN
        self.classID.objectName = "Polymer"

#######################################################
#                 GETTER METHODS                      #
#######################################################
    def getPos(self):
        temp = []
        for i in range(len(self.positions)):
            temp.append(self.positions[i].copy())
        return temp
    def getBonds(self,n=0,i=0):
        #print(offset)

        temp = []
        offset = n*i + i
        for i in range(len(self.bonds)):
            temp.append([])
            for j in range(len(self.bonds[i])):
                #print(temp[i][j])
                temp[i].append(self.bonds[i][j] + offset)
        #print(temp)
        return temp
    def getIDs(self):
        return self.ids.copy()


#######################################################


    def createPolymerChain(self,N=10,rez=1,deg = [0,0], origin=[0,0,0],func=None,rnd = 100):
        if func is None:
            deg[0]=deg[0]*math.pi/180
            deg[1]=deg[1]*math.pi/180
            for i in range(N):
                x = round(math.cos(deg[0])*math.sin(deg[1])*i*rez,rnd) + origin[0]
                y = round(math.sin(deg[0])*math.sin(deg[1])*i*rez,rnd) + origin[1]
                z = round(math.cos(deg[0])*i*rez,rnd) + origin[2]
                self.positions.append([x,y,z])

        elif func is not None:
            for i in range(N):
                if func(i,N) is not None:
                    self.positions.append(add(func(i,N),origin))


    def move(self,dest,origin='G'):
        if origin is 'G':

            for i in range(len(self.positions)):
                    self.positions[i][0] =dest[0]
                    self.positions[i][1] =dest[1]
                    self.positions[i][2] =dest[2]


        if origin is 'L':

            for i in range(len(self.positions)):
                    self.positions[i][0] +=dest[0]
                    self.positions[i][1] +=dest[1]
                    self.positions[i][2] +=dest[2]



    def defineBonds(self,func=None,custom=None,offset=0):
        N = len(self.positions)

        if func != None:
            for i in range(N):
                if func(i,N) is not None:
                    self.bonds.append(func(i,N))

        ########################
        elif custom != None:
            self.bonds = custom

        ########################
        else:
            for i in range(len(self.positions)-1):
                self.bonds.append([i,i+1])
        return None




    def defineParticleTypes(self,n=1,func=None,custom=None):
        N = len(self.positions)

        if func != None:
            for i in range(N):
                if func(i,N) is not None:
                    self.ids.append(func(i,N))

        ########################
        elif custom != None:
            self.ids = custom
            return None

        ########################
        else:
            t = 0
            for i in range(N):
                if (i % n ==0):
                    t=0
                ids.append(t)
                t = t+1
        return None
