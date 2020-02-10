import Simulations

parameters = Simulations.PolymerSimulationParameters()

parameters.setLength(10)
parameters.setNumberChains(1)
parameters.setPairRadius(0.3)
parameters.setPairPotentialStrength(10**3)
parameters.setPairRadiusEqualibrium(0.3)
parameters.setPairMaximumRadius(0.1)
parameters.setForcePull(10)
parameters.setPeriodicAmplitude(10)
parameters.setGamma(0.5)
parameters.setKBT(1)
parameters.setTimeStep(0.000001)
parameters.setProbePeriod(10000)
parameters.setRunLength(1000000)

sim = Simulations.PolymerSimulation(parameter=parameters)

filelocation = sim.run()
print(filelocation)
renderer = Simulations.DataVisualizer("Simulation_2020-02-07_22-22-28/Simulation_2020-02-07_22-22-28",interval=0)
