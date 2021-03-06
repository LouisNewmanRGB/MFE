from SimulationPython.Simulation.Particle3D import Particle3D
from SimulationPython.Simulation.Environment import Environment
from SimulationPython.Simulation.Simulation import Simulation
from SimulationPython.Simulation.Sphere import Sphere
from SimulationPython.Simulation.Ellipsoid import Ellipsoid

nStep = 10
nPart = 2
timeStep = 0.1 #ms
D = 2 #um2/ms
l = (6*D*timeStep)**0.5
envSize = 10*l
T2 = 1 #irrelevant
#perm = 0.034/1000 #mm/ms
probInOut = 0.5

envi = Environment(T2, D, envSize, envSize, envSize)
part = [Particle3D(0.,0.,0.) for i in range(nPart)]
compartments = [Sphere(0,0,0,T2,2*D, probInOut, 2*l), Ellipsoid(0, 0, 0, T2, 2*D, probInOut, 2*l, 3*l, 1*l)]
#compartments = [Ellipsoid(0, 0, 0, T2, 2*D, perm, np.array([2*l, 2*l, 0]), np.array([l, -l, 0]), np.array([0, 0, 3*l]))]
sim = Simulation(nStep, timeStep, part, envi, compartments)
sim.run(seed=None, calcData = True)

sim.plot(False)
sim.plot(True)
print(sim.getStepLengths())
