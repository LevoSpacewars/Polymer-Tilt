import Simulations
import math
# parameters = Simulations.PolymerSimulationParameters()
# dt = 0.000001
# runl = 10000000
# vAverage = 0.05
# vs = 1000000
# A=1
# print("\n\n\n\n\n")
# parameters.setSheerForceRange(0,1.5)
# parameters.setDf(1)
# parameters.setLength(100)
# parameters.setNumberChains(10)
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
# sim.init(parameter=parameters,initializer='--mode=gpu')
#
# filelocation = sim.run()
# print(filelocation)

renderer = Simulations.DataVisualizer(basedirectory="Simulation_2020-05-17_23-48-14/",interval=0)
renderer.init(plotTilt=False,plotProbabilityMap = True,animatePolymers=True )

## TODO: implement stderr of mean to see if mean is still moving
