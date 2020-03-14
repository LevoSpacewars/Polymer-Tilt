import Simulations
import math
parameters = Simulations.PolymerSimulationParameters()
dt = 0.0000001
runl = 1000000
vAverage = 0.05
vs = 90000
A=1
print("\n\n\n\n\n")
parameters.setLength(10)
parameters.setNumberChains(2)
parameters.setPairRadius(.3)
parameters.setPairPotentialStrength(10e4)
parameters.setPairRadiusEqualibrium(.3)
parameters.setPairMaximumRadius(0.03)
parameters.setForcePull(10)
parameters.setPeriodicAmplitude(A)
parameters.setGamma(0.5)
parameters.setKBT(1)
parameters.setTimeStep(dt)
parameters.setProbePeriod(vs)
parameters.setRunLength(runl)

sim = Simulations.PolymerSimulation(parameter=parameters)

filelocation = sim.run(forceRange=[10,11,1])
print(filelocation)
renderer = Simulations.DataVisualizer(basedirectory=filelocation,interval=0)
