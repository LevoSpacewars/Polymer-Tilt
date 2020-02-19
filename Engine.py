import Simulations


parameters = Simulations.PolymerSimulationParameters()
dt = 0.00001
runl = 10000000
vAverage = 0.05
vs = vAverage/dt
print(runl/vs)
parameters.setLength(30)
parameters.setNumberChains(3)
parameters.setPairRadius(0.3)
parameters.setPairPotentialStrength(10**3)
parameters.setPairRadiusEqualibrium(0.3)
parameters.setPairMaximumRadius(0.2)
parameters.setForcePull(10)
parameters.setPeriodicAmplitude(10)
parameters.setGamma(0.5)
parameters.setKBT(1)
parameters.setTimeStep(dt)
parameters.setProbePeriod(vs)
parameters.setRunLength(runl)

sim = Simulations.PolymerSimulation(parameter=parameters)

filelocation = sim.run(forceRange=[0,10,3])
print(filelocation)
renderer = Simulations.DataVisualizer(basedirectory=filelocation)
