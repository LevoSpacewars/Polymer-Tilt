import Simulations

parameters = Simulations.PolymerSimulationParameters()

parameters.setLength(10)
parameters.setNumberChains(10)
parameters.setPairRadius(0.3)
parameters.setPairPotentialStrength(10**3)
parameters.setPairRadiusEqualibrium(0.3)
parameters.setPairMaximumRadius(0.1)
parameters.setForcePull(10)
parameters.setPeriodicAmplitude(10)
parameters.setGamma(0.5)
parameters.setKBT(1)
parameters.setTimeStep(0.000001)
parameters.setProbePeriod(1000)
parameters.setRunLength(1000000)

sim = Simulations.PolymerSimulation(parameter=parameters)

filelocation = sim.run()
renderer = Simulations.DataVisualizer(filelocation)
