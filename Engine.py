import Simulations
import math
parameters = Simulations.PolymerSimulationParameters()
dt = 0.000001
runl = 100
vs = 50
A=1
print("\n\n\n\n\n")
parameters.setLength(100)
parameters.setNumberChains(8)
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

filelocation = sim.run(forceRange=[0,3,20])
print(filelocation)

renderer = Simulations.DataVisualizer(basedirectory=filelocation,interval=0)
renderer.init()
a = Simulations.GlobalDataAnalyzer(filelocation)

a.plotTiltbyForce(20,1)
