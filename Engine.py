import Simulations
import math
parameters = Simulations.PolymerSimulationParameters()
dt = 0.001
runl =  2*10**6
vs =    1*10**0
A= 0.3/0.1 * 1
print("\n\n\n\n\n")

parameters.setSheerForceRange(0,2)
parameters.setDf(15)
parameters.setLength(100)
parameters.setNumberChains(10)
parameters.setPairRadius(0.1)
parameters.setPairPotentialStrength(10e3)
parameters.setPairRadiusEqualibrium(0.1)
parameters.setPairMaximumRadius(0.1)
parameters.setForcePull(20)
parameters.setPeriodicAmplitude(A)
parameters.setGamma(0.5)
parameters.setKBT(1)
parameters.setTimeStep(dt)
parameters.setProbePeriod(vs)
parameters.setRunLength(runl)
parameters.setIntegrator("legavin")
parameters.setRunDirection("forward")

sim = Simulations.PolymerSimulation()
sim.init(parameter=parameters,initializer='--mode=cpu')

filelocation = sim.run(server = True)
print(filelocation)

#renderer = Simulations.DataVisualizer(basedirectory=filelocation,interval=0.5)
#renderer.init(plotTilt=False,plotProbabilityMap = True,animatePolymers=False,plotPolymerProfiles=False)
