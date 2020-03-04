import Simulations

parameters = Simulations.PolymerSimulationParameters()
dt = 0.0000001
runl = 100000000
vAverage = 0.05
vs = 1000000
print("\n\n\n\n\n")
parameters.setLength(100)
parameters.setNumberChains(8)
parameters.setPairRadius(.3)
parameters.setPairPotentialStrength(10e4)
parameters.setPairRadiusEqualibrium(.3)
parameters.setPairMaximumRadius(0.03)
parameters.setForcePull(10)
parameters.setPeriodicAmplitude(10)
parameters.setGamma(0.5)
parameters.setKBT(1)
parameters.setTimeStep(dt)
parameters.setProbePeriod(vs)
parameters.setRunLength(runl)

sim = Simulations.PolymerSimulation(parameter=parameters)

filelocation = sim.run(forceRange=[0,0,1])
print(filelocation)
renderer = Simulations.DataVisualizer(basedirectory=filelocation,interval=0)
