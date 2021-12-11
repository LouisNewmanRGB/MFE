from SimulationPython.Simulation.Particle3D import Particle3D
from SimulationPython.Simulation.Environment import Environment
from SimulationPython.Simulation.Simulation import Simulation
from SimulationPython.Simulation.Sphere import Sphere

nStep = 8
nPart = 1
timeStep = 1 #ms
#timeStep = 1
D = 2 #um2/ms
l = (6*D*timeStep)**0.5
radius = 8 #um
envSize = 3*radius
T2 = 1 #irrelevant
#perm = 0.034/1000 #mm/ms
permeability = 0.5

envi = Environment(T2, D, envSize, envSize, envSize)
#part = [Particle3D(*Util.getRandomDirection()*Util.getRandomQuadratic(radius)) for i in range(1)]
part = [Particle3D(0, 0, 0) for i in range(nPart)]
compartments = [Sphere(0,0,0,T2,2*D, permeability, radius)]
sim = Simulation(nStep, timeStep, part, envi, compartments)

sim.run(seed=None, calcData = True)

sim.plot(False)
sim.plot(True)
print(sim.getStepLengths())
#print(sim.getPositions())
#plt.hist([np.linalg.norm(pos) for pos in sim.getPositions()], bins=20)
#plt.show()
