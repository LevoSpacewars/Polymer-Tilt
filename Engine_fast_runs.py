import Simulations
import math
import sys
parameters = Simulations.PolymerSimulationParameters()

print("\n\n\n\n\n")
#sys.args = [length, kbt, amplitude, tension, #polymers, runtime, samplerate, sheer, ID]



length = int(sys.argv[1])
kbt = float(sys.argv[2])
phi0 = float(sys.argv[3])
ten = float(sys.argv[4])
chainnum = int(sys.argv[5])
time_total = int(sys.argv[6])
sample_rate = int(sys.argv[7])
sheer = float(sys.argv[8])
name = sys.argv[9]

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

sim = Simulations.PolymerSimulation()
sim.init(parameter=parameters,initializer='--mode=gpu')



filelocation = sim.probe(name,sheer,"",server = False) ##S
print(filelocation)

#renderer = Simulations.DataVisualizer(basedirectory=filelocation,interval=0.75)
#renderer.init(plotTilt=False,plotProbabilityMap = True,animatePolymers=False,plotPolymerProfiles=False)
