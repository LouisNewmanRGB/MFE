import numpy as np
import matplotlib.pyplot as plt

from Particle3D import Particle3D
from Environment import Environment
from Simulation import Simulation
from Sphere import Sphere

nStep = 8
timeStep = 10/nStep #ms
D = 2 #um2/ms
l = (6*D*timeStep)**0.5
radius = 8 #um
envSize = 3*radius
T2 = 1 #irrelevant
#perm = 0.034/1000 #mm/ms
probInOut = 0.5

envi = Environment(T2, D, envSize, envSize, envSize)
part = [Particle3D(0.,0.,0.) for i in range(2)]
compartments = [Sphere(0,0,0,T2,2*D, probInOut, radius)]
sim = Simulation(nStep, timeStep, part, envi, compartments)

sim.run(seed=None, calcData = True)

sim.plot(False)
sim.plot(True)
print(sim.getStepLengths())
