import numpy as np
import matplotlib.pyplot as plt

from SimulationCpp import Simulation
from SimulationPython.Simulation.Util import Util
from SimulationCppResults import SimulationCppResults

nStep = 8
nPart = 5*int(1e4)
diffusionTime = 50
delta = 5
Delta = diffusionTime
TE = 2*diffusionTime
T2 = np.inf
D = 2

radius = 8
envSize = 5*radius

sim = Simulation(nStep, 1.1*(delta+Delta)/nStep, T2, D, envSize, envSize, envSize)
#sim.addSphere(0, 0, 0, T2, D, 0, radius)
#sim.addEllipsoid(0, 0, 0, T2, D, 0, radius*np.array([1, 0, 0]), radius*np.array([0, 1, 0]), radius*np.array([0, 0, 1]))
sim.addEllipsoid(0, 0, 0, T2, D, 0, radius, radius, radius)
sim.createStartingPositions(nPart, True)
sim.createSequencePWG(np.array([delta, Delta, Delta+delta, TE]), np.array([[-1, 0, 0], [0, 0, 0], [1, 0, 0], [0, 0, 0]]))
#sim.createSequenceSGP()
sim.run()
dataVector = sim.getResults()
results = SimulationCppResults(dataVector, nPart)

if nPart <= 500:
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    #ax.scatter(results.startingPositions[:, 0], results.startingPositions[:, 1], results.startingPositions[:, 2])
    ax.scatter(results.getEndPositions()[:, 0], results.getEndPositions()[:, 1], results.getEndPositions()[:, 2])
    #ax.scatter(results.endPositionsVoxel[:, 0], results.endPositionsVoxel[:, 1], results.endPositionsVoxel[:, 2])
    plt.show()

qPoints = np.linspace(0, 0.4, 101)[1:]
qVectors = np.array([[q, 0, 0] for q in qPoints])
#simulatedSignal, stdsSample, stdsPopulation = results.getSGPSignal(qVectors, includeStd=True)
simulatedSignal = np.empty(len(qPoints))
stdsSample = np.empty(len(qPoints))
stdsPopulation = np.empty(len(qPoints))
for q in range(len(qPoints)):
    simulatedSignal[q], stdsSample[q], stdsPopulation[q] = results.getFiniteGradientSignalGammaG(qPoints[q]/delta, includeStd=True)
trueSignalPoints = Util.getSignal_sphere_fin(radius, D, diffusionTime, 20, 20, 10)(qPoints)

plt.plot(qPoints, trueSignalPoints, color="r")#colors[t], marker=".")
plt.errorbar(qPoints, simulatedSignal, stdsSample, fmt=".", color="g")
plt.legend(["Theoretical signal", "Simulated signal", "Simulated signal (real part)"])
plt.xlabel("q=gamma*G*delta [um-1]")
plt.ylabel("Signal attenuation")
plt.title("graphTitle")
plt.yscale("log")
plt.grid()
plt.show()
