import numpy as np
import matplotlib.pyplot as plt

from Particle3D import Particle3D
from Environment import Environment
from Simulation import Simulation

nPart = 2
nStep = 8
timeStep = 10/nStep #ms
timeStep = 1
D = 2 #um2/ms
l = (6*D*timeStep)**0.5
envSize = 5*l
T2 = 1 #irrelevant

envi = Environment(T2, D, envSize, envSize, envSize)
part = [Particle3D(0, 0, 0) for i in range(nPart)]
sim = Simulation(nStep, timeStep, part, envi)
sim.run(seed=None, calcData = True)

sim.plot(False)
sim.plot(True)
