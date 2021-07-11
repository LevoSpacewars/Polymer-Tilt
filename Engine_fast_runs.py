from matplotlib.pyplot import legend
import Simulations
import math
import random as r
import sys
from disorder import get_amount_nodes
parameters = Simulations.PolymerSimulationParameters()

print("\n\n\n\n\n")
inputs = {"length": None, "kbt": None, "amplitude": None, "tension": None, "npolymers": None, "runtime": None, "samplerate": None, "sheer": None, "endsheer": None, "random_seed": None, "disorder_ratio": None, "ID": None}

# TODO: need to implement the disorder level to set the contraints along the min and max amplitude


print(sys.argv)
length = int(sys.argv[1])
kbt = float(sys.argv[2])
phi0 = float(sys.argv[3])
ten = float(sys.argv[4])
chainnum = int(sys.argv[5])
time_total = int(sys.argv[6])
sample_rate = int(sys.argv[7])
sheer = float(sys.argv[8])
endsheer = float(sys.argv[9])
id = sys.argv[10]
random_seed = 0
disorder_level = 0
if (len(sys.argv) > 11):
    random_seed = int(sys.argv[11])
    amount_nodes = get_amount_nodes (float(sys.argv[12]))





dt = 0.001
runl = time_total
vs =    sample_rate
A= phi0


parameters.setSheerForceRange(sheer,0)
parameters.setDf(1)
parameters.setLength(length)
parameters.setNumberChains(chainnum)
parameters.setPairRadius(0.1)
parameters.setPairPotentialStrength(10e3)
parameters.setPairRadiusEqualibrium(0.1)
parameters.setPairMaximumRadius(0.1)
parameters.setForcePull(ten)
parameters.setPeriodicAmplitude(phi0)
parameters.setGamma(0.5)
parameters.setKBT(kbt)
parameters.setTimeStep(dt)
parameters.setProbePeriod(vs)
parameters.setRunLength(runl)
parameters.setIntegrator("legavin")
parameters.setRunDirection("forward")

amplitude_range = (0.005* phi0, 0.05 * phi0)
sim = Simulations.PolymerSimulation()
sim.set_disorder(random_seed,amplitude_range,amount_nodes,chainnum)
print(amount_nodes)
sim.init(parameter=parameters,initializer='--mode=gpu')


filelocation = sim.probe(id,sheer,"",server = True) ##S
print(filelocation)

#renderer = Simulations.DataVisualizer(basedirectory=filelocation,interval=0.75)
#renderer.init(plotTilt=False,plotProbabilityMap = True,animatePolymers=False,plotPolymerProfiles=False)
