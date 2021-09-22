import numpy as np

from Sphere import Sphere
from Environment import Environment
from Particle3D import Particle3D
from Simulation import Simulation
from Util import Util

nSpheres = 5000
radius = 8#um
voxelDim = 500#um
nStep = 8
nPart = int(1e5)
D = 2e-3*1e6*1e-3 #um2/ms
probInOut = 0.25
T2 = 1 #not useful

diffusionTimes = np.array([1, 20, 100])#ms
timeSteps = diffusionTimes/nStep
for timeStep in timeSteps:
    l = (6*D*timeStep)**0.5
    print("Time step = {timeStep}ms\nl = {l}um\n".format(timeStep=timeStep, l=l))

timeStep = timeSteps[0]
print(timeStep)

nSpheresLine = int((5000)**(1/3))
#print(nSpheresLine)
sphereDistance = voxelDim/nSpheresLine - 2*radius
#print(sphereDistance)
spacing = voxelDim / (2*nSpheresLine)
#print(spacing)
#print(sphereDistance/2 + radius)
positions = np.linspace(spacing, voxelDim - spacing, nSpheresLine)
#print(positions)
#print([(2*i + 1)*spacing for i in range(nSpheresLine)])


compartments = []
for x in positions:
    for y in positions:
        for z in positions:
            compartments.append(Sphere(x, y, z, T2, D, probInOut, radius))
env = Environment(0, 0, 0, T2, D, voxelDim, voxelDim, voxelDim)
particles = [Particle3D(Util.getRandomU(voxelDim),Util.getRandomU(voxelDim),Util.getRandomU(voxelDim)) for i in range(nPart)]
sim = Simulation(nStep, timeStep, particles, env, compartments)
sim.run(seed=1, calcData=False)
print("Done")
