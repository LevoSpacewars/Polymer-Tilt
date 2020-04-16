import Simulations
#import math
#parameters = Simulations.PolymerSimulationParameters()
#dt = 0.000001
#runl = 100000
#vAverage = 0.05
#vs = 100
#A=1
#print("\n\n\n\n\n")
#parameters.setSheerForceRange(0,2)
#parameters.setDf(10)
#parameters.setLength(100)
#parameters.setNumberChains(10)
#parameters.setPairRadius(.3)
#parameters.setPairPotentialStrength(10e4)
#parameters.setPairRadiusEqualibrium(.3)
#parameters.setPairMaximumRadius(0.03)
#parameters.setForcePull(20)
#parameters.setPeriodicAmplitude(A)
#parameters.setGamma(0.5)
#parameters.setKBT(1)
#parameters.setTimeStep(dt)
#parameters.setProbePeriod(vs)
#parameters.setRunLength(runl)

#sim = Simulations.PolymerSimulation()
#sim.init(parameter=parameters,initializer='--mode=gpu')

#filelocation = sim.runTest()
#print(filelocation)
#
renderer = Simulations.DataVisualizer(basedirectory="Simulation_2020-04-08_13-25-13/",interval=0)
renderer.init()
