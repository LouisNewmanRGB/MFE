import numpy as np
import matplotlib.pyplot as plt

from Particle3D import Particle3D
from Environment import Environment
from Simulation import Simulation

# nStep = 10
# timeStep = 1
# D = 1
# part = [Particle3D(0.,0.,0.) for i in range(500)]
# envi = Environment(0, 0, 0, 1, D, 1, 1, 1)
#
# sim = Simulation(nStep, timeStep, part, envi, [])
#
# sim.run(seed=None, calcData = False)
#
# plt.hist(sim.getDistances())
# print(np.average(np.power(sim.getDistances(), 2))) #mean squared distance
# print(6*D*nStep*timeStep) #theoretical mean squared distance
# plt.show()

nStep = 10
timeStep = 1
D = 1
part = [Particle3D(0.,0.,0.) for i in range(2)]
envi = Environment(0, 0, 0, 1, D, 5, 5, 5)

sim = Simulation(10, 1, part, envi, [])

sim.run(seed=None, calcData = True)

sim.plot(False)
sim.plot(True)