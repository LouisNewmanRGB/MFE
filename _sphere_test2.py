import numpy as np
import matplotlib.pyplot as plt

from Particle3D import Particle3D
from Environment import Environment
from Simulation import Simulation
from Sphere import Sphere

nStep = 10
timeStep = 1
D = 1
part = [Particle3D(0.,0.,0.) for i in range(2)]
envi = Environment(0, 0, 0, 1, D, 5, 5, 5)
compartments = [Sphere(0,0,0,2,0.5*D,2)]
sim = Simulation(10, 1, part, envi, compartments)

sim.run(seed=1, calcData = True)

sim.plot(False)
sim.plot(True)
print(sim.getStepLengths())
