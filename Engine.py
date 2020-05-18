import Simulations
import math
parameters = Simulations.PolymerSimulationParameters()
dt = 0.000001
runl = 100000000
vAverage = 0.05
vs = 1000000
A=1
print("\n\n\n\n\n")
parameters.setSheerForceRange(0,1.5)
parameters.setDf(10)
parameters.setLength(200)
parameters.setNumberChains(10)
parameters.setPairRadius(.3)
parameters.setPairPotentialStrength(10e4)
parameters.setPairRadiusEqualibrium(.3)
parameters.setPairMaximumRadius(0.03)
parameters.setForcePull(20)
parameters.setPeriodicAmplitude(A)
parameters.setGamma(0.5)
parameters.setKBT(1)
parameters.setTimeStep(dt)
parameters.setProbePeriod(vs)
parameters.setRunLength(runl)

sim = Simulations.PolymerSimulation()
sim.init(parameter=parameters,initializer='--mode=gpu')

filelocation = sim.run()
print(filelocation)

renderer = Simulations.DataVisualizer(basedirectory=filelocation,interval=0)
renderer.init(plotTilt=True,plotProbabilityMap = True )
