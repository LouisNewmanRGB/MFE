import numpy as np

from SimulationPython.Simulation.Sphere import Sphere
from SimulationPython.Simulation.Environment import Environment
from SimulationPython.Simulation.Particle3D import Particle3D
from SimulationPython.Simulation.Util import Util
from SimulationPython.Simulation.SimulationParallel import SimulationParallel

if __name__ == "__main__":

    nSpheres = 5000
    radius = 8#um
    voxelDim = 500#um
    nStep = 8
    nPart = int(1e5)
    D = 2 #um2/ms
    probInOut = 0.25
    T2 = 1 #not useful

    diffusionTimes = np.array([1, 20, 100])#ms
    timeSteps = diffusionTimes/nStep
    for timeStep in timeSteps:
        l = (6*D*timeStep)**0.5
        print("Time step = {timeStep}ms\nl = {l}um\n".format(timeStep=timeStep, l=l))

    timeStep = timeSteps[1]
    print("Chosen timestep: {timeStep}ms".format(timeStep=timeStep))

    nSpheresLine = int((5000)**(1/3))
    #sphereDistance = voxelDim/nSpheresLine - 2*radius
    spacing = voxelDim / (2*nSpheresLine)
    positions = np.linspace(spacing, voxelDim - spacing, nSpheresLine)

    def getCompartments():
        compartments = []
        for x in positions:
            for y in positions:
                for z in positions:
                    compartments.append(Sphere(x, y, z, T2, D, probInOut, radius))
        return compartments
    env = Environment(T2, D, voxelDim, voxelDim, voxelDim)
    particles = [Particle3D(Util.getRandomU(voxelDim),Util.getRandomU(voxelDim),Util.getRandomU(voxelDim)) for i in range(nPart)]
    sim = SimulationParallel(nStep, timeStep, particles, env, getCompartments)
    sim.run(seed=1, calcData=False, partPrintNumber=10)
    print("Done")
