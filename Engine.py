import Simulations


parameters = Simulations.PolymerSimulationParameters()
dt = 0.0000001
runl = 10000000
vAverage = 0.05
vs = 100000
print("\n\n\n\n\n")
parameters.setLength(10)
parameters.setNumberChains(1)
parameters.setPairRadius(.3)
parameters.setPairPotentialStrength(10**3)
parameters.setPairRadiusEqualibrium(.3)
parameters.setPairMaximumRadius(.1)
parameters.setForcePull(10)
parameters.setPeriodicAmplitude(1)
parameters.setGamma(0.5)
parameters.setKBT(1)
parameters.setTimeStep(dt)
parameters.setProbePeriod(vs)
parameters.setRunLength(runl)

sim = Simulations.PolymerSimulation(parameter=parameters)

filelocation = sim.run(forceRange=[0,0,1])
print(filelocation)
renderer = Simulations.DataVisualizer(basedirectory=filelocation)
