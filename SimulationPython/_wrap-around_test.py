from SimulationPython.Simulation.Particle3D import Particle3D
from SimulationPython.Simulation.Environment import Environment
from SimulationPython.Simulation.Simulation import Simulation

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

print(sim.getStepLengths())
sim.plot(False)
sim.plot(True)
