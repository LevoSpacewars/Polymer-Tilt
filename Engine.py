import Simulations


parameters = Simulations.PolymerSimulationParameters()
dt = 0.000001
runl = 10000000
vAverage = 0.05
vs = vAverage/dt
print(runl/vs)
print("\n\n\n\n\n")
parameters.setLength(50)
parameters.setNumberChains(10)
parameters.setPairRadius(0.3)
parameters.setPairPotentialStrength(10**3)
parameters.setPairRadiusEqualibrium(0.3)
parameters.setPairMaximumRadius(0.03)
parameters.setForcePull(100)
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
