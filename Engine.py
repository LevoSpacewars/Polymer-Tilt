import Simulations

#
# parameters = Simulations.PolymerSimulationParameters()
# dt = 0.00001
# runl = 10000000
# vAverage = 0.5
# vs = 0.5/dt
# print(runl/vs)
# parameters.setLength(50)
# parameters.setNumberChains(50)
# parameters.setPairRadius(0.3)
# parameters.setPairPotentialStrength(10**3)
# parameters.setPairRadiusEqualibrium(0.3)
# parameters.setPairMaximumRadius(0.1)
# parameters.setForcePull(10)
# parameters.setPeriodicAmplitude(10)
# parameters.setGamma(0.5)
# parameters.setKBT(1)
# parameters.setTimeStep(0.000001)
# parameters.setProbePeriod(0.5/dt)
# parameters.setRunLength(runl)
#
# sim = Simulations.PolymerSimulation(parameter=parameters)
#
# filelocation = sim.run(forceRange=[0,10,10])
# print(filelocation)
renderer = Simulations.DataVisualizer(directory="",gsdname="polymer_0.0",parametername="polymer_0.0_simulation_parameters.txt",interval=0)
