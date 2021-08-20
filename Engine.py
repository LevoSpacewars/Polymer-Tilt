import Simulations
import math
from random import randint
from disorder import get_amplitude_mod
parameters = Simulations.PolymerSimulationParameters()
dt = 0.001
runl = (1*10**6)
vs =    1*10**3
A= -0.3/0.1 * 0.0001
print("\n\n\n\n\n")

parameters.setSheerForceRange(0,4)
parameters.setDf(10)
parameters.setLength(50)
parameters.setNumberChains(1)
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
parameters.setRunDirection("forward")
parameters.setDisorder(0)
parameters.setCommensurate(True)
mod = get_amplitude_mod(0.600,60,abs(-3))

amplitude_range = (0.005, 0.005 * mod)
print(amplitude_range)

sim = Simulations.PolymerSimulation()
sim.init(parameter=parameters,initializer='--mode=cpu')
sim.set_disorder(randint(0,199999),amplitude_range,60,10,0,3)
filelocation = sim.run(name="singlePolymert")
print(filelocation)

#renderer = Simulations.DataVisualizer(basedirectory=filelocation,interval=0.75)
#renderer.init(plotTilt=False,plotProbabilityMap = True,animatePolymers=False,plotPolymerProfiles=False)
