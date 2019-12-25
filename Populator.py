import hoomd
import hoomd.md
import Util
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
        return self.positions
    def getBonds(self):
        return self.bonds
    def getIDs(self):
        return self.ids


#######################################################


    def createPolymerChain(self,N=10,di = [0,1,0], origin=[0,0,0],func=None):
        if func is None:
            for i in range(N):
                x = di[0]*i + origin[0]
                y = di[1]*i + origin[1]
                z = di[2]*i + origin[2]
                self.positions.append([x,y,z])

        elif func is not None:
            for i in range(N):
                if func(i,N) is not None:
                    self.positions.append(add(func(i,N),origin))




    def defineBonds(self,func=None,custom=None):
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
            for i in range(1,len(positions)):
                self.bonds.append([i-1,i])
        return None




    def defineParticleTypes(n,func=None):
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
            t = 0
            for i in range(N):
                if i % n ==0:
                    t=0
                ids.append(t++)
        return None
class Populator:
