import numpy as np
import matplotlib.pyplot as plt

from SimulationCpp import Simulation
from SimulationPython.Simulation.Util import Util
from SimulationCppResults import SimulationCppResults

nStep = 8
nPart = int(1e6)
diffusionTime = 50
delta = 5
Delta = diffusionTime
TE = 2*diffusionTime
T2 = np.inf
D = 2

radius = 8
envSize = 5*radius

sim = Simulation(nStep, Delta/nStep, T2, D, envSize, envSize, envSize)

sim.addCylinderBasic(radius, radius/2, T2, D, 0, radius)
sim.createStartingPositions(nPart, True)
#sim.createSequencePWG(np.array([delta, Delta, Delta+delta, TE]), np.array([[-1, 0, 0], [0, 0, 0], [1, 0, 0], [0, 0, 0]]))
sim.createSequenceSGP()
sim.run()
dataVector = sim.getResults()
results = SimulationCppResults(dataVector, nPart)

if nPart <= 1000:
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(results.startingPositions[:, 0], results.startingPositions[:, 1], results.startingPositions[:, 2])
    #ax.scatter(results.getEndPositions()[:, 0], results.getEndPositions()[:, 1], results.getEndPositions()[:, 2])
    #ax.scatter(results.endPositionsVoxel[:, 0], results.endPositionsVoxel[:, 1], results.endPositionsVoxel[:, 2])
    plt.show()

qValues = np.linspace(0, 0.8, 101)[1:]
qVectors = np.array([[q, 0, 0] for q in qValues])
signal, stdsSample, stdsPopulation = results.getSGPSignal(qVectors, includeStd=True)
trueSignalFunc = Util.getSignal_cylinder_fin(radius, D, Delta, 20, 20, 5)
trueSignalPoints = trueSignalFunc(qValues)

plt.yscale("log")
plt.errorbar(qValues, signal, stdsSample)
plt.plot(qValues, trueSignalPoints, color="r")
plt.show()
