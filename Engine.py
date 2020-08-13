import Simulations
import math
parameters = Simulations.PolymerSimulationParameters()
dt = 0.0001
runl =  1*10**7
vs =    1*10**4
A= 0.3/0.1 * 1
print("\n\n\n\n\n")
## TODO:  add move functionality for the 0frame saves
# TODO: Change how the formatter formats, So that the force for a value is coupled. This way I can merge other graphs into a single output file.
parameters.setSheerForceRange(0,2)
parameters.setDf(30)
parameters.setLength(200)
parameters.setNumberChains(10)
parameters.setPairRadius(0.1)
parameters.setPairPotentialStrength(10e3)
parameters.setPairRadiusEqualibrium(0.1)
parameters.setPairMaximumRadius(0.1)
parameters.setForcePull(20)
parameters.setPeriodicAmplitude(A)
parameters.setGamma(0.5)
parameters.setKBT(1)
parameters.setTimeStep(dt)
parameters.setProbePeriod(vs)
parameters.setRunLength(runl)
parameters.setIntegrator("legavin")

sim = Simulations.PolymerSimulation()
sim.init(parameter=parameters,initializer='--mode=gpu')

filelocation = sim.run()
print(filelocation)

renderer = Simulations.DataVisualizer(basedirectory=filelocation,interval=0.9)
renderer.init(plotTilt=True,plotProbabilityMap = True,animatePolymers=False,plotPolymerProfiles=False)
