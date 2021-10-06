import numpy as np
import matplotlib.pyplot as plt

from Particle3D import Particle3D
from Environment import Environment
from Simulation import Simulation
from Sphere import Sphere
from Util import Util

nStep = 2
timeStep = 100/nStep #ms
#timeStep = 1
D = 2 #um2/ms
l = (6*D*timeStep)**0.5
radius = 8 #um
envSize = 3*radius
T2 = 1 #irrelevant
#perm = 0.034/1000 #mm/ms
probInOut = 0

envi = Environment(T2, D, envSize, envSize, envSize)
#part = [Particle3D(*Util.getRandomDirection()*Util.getRandomQuadratic(radius)) for i in range(1)]
part = [Particle3D(0, 0, 0) for i in range(1)]
compartments = [Sphere(0,0,0,T2,2*D, probInOut, radius)]
sim = Simulation(nStep, timeStep, part, envi, compartments)

sim.run(seed=None, calcData = True)

sim.plot(False)
sim.plot(True)
print(sim.getStepLengths())
#print(sim.getPositions())
#plt.hist([np.linalg.norm(pos) for pos in sim.getPositions()], bins=20)
#plt.show()
