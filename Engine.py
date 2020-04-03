import Simulations
import math
# parameters = Simulations.PolymerSimulationParameters()
# dt = 0.0000001
# runl = 100000000
# vAverage = 0.05
# vs = 1000000
# A=1
# print("\n\n\n\n\n")
# parameters.setLength(100)
# parameters.setNumberChains(8)
# parameters.setPairRadius(.3)
# parameters.setPairPotentialStrength(10e4)
# parameters.setPairRadiusEqualibrium(.3)
# parameters.setPairMaximumRadius(0.03)
# parameters.setForcePull(20)
# parameters.setPeriodicAmplitude(A)
# parameters.setGamma(0.5)
# parameters.setKBT(1)
# parameters.setTimeStep(dt)
# parameters.setProbePeriod(vs)
# parameters.setRunLength(runl)
#
# sim = Simulations.PolymerSimulation()
#
# sim.init(parameter=parameters,initializer='--mode=gpu')
#
# filelocation = sim.run(forceRange=[0,5,10])
# print(filelocation)

#renderer = Simulations.DataVisualizer(basedirectory="Simulation_2020-04-01_22-01-19/",interval=0.5)
#renderer.init()
filelocation="Simulation_2020-04-01_22-01-19/"
a = Simulations.GlobalDataAnalyzer(filelocation)

a.plotTiltbyForce(20,1)
