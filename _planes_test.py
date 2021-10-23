import numpy as np
import matplotlib.pyplot as plt

from Particle3D import Particle3D
from Environment import Environment
from Simulation import Simulation
from Planes import Planes
from Util import Util

nPart = 2
nStep = 3
#timeStep = 10/nStep #ms
timeStep = 10
D = 2 #um2/ms
l = (6*D*timeStep)**0.5
envSize = 5*l
spacing = 2*l
T2 = 1 #irrelevant

envi = Environment(T2, D, envSize, envSize, envSize)
part = [Particle3D(Util.getRandomU(spacing), Util.getRandomU(envSize), Util.getRandomU(envSize)) for i in range(nPart)]
sim = Simulation(nStep, timeStep, part, envi, [Planes(T2, D, spacing)])
sim.run(seed=None, calcData = True)

print(sim.getStepLengths())
sim.plot(False)
sim.plot(True)
