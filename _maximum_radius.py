import numpy as np
import matplotlib.pyplot as plt

from SimulationPython.Simulation.Util import Util

D = 2
diffusionTime = 400 #50
spacing = 40 #um
qPoints = np.linspace(0, 0.2, 101)[1:]

"""
#data = np.load("./RandomWalkSimulator/results/npyTest.npy")
data = np.load("./RandomWalkSimulator/results/planes_max_radius.npy")
print(data.shape)
print(data[0:7, :])

results = SimulationResults2(data)
if data.shape[0] < 1000:
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    #ax.scatter(results.startingPositions[:, 0], results.startingPositions[:, 1], results.startingPositions[:, 2])
    #ax.scatter(results.getEndPositions()[:, 0], results.getEndPositions()[:, 1], results.getEndPositions()[:, 2])
    ax.scatter(results.endPositionsVoxel[:, 0], results.endPositionsVoxel[:, 1], results.endPositionsVoxel[:, 2])
    plt.show()

qVectors = np.array([[q, 0, 0] for q in qPoints])
simulatedSignal, stdsSample, stdsPopulation = results.getSGPSignal(qVectors, includeStd=True)

trueSignalPoints = Util.getSignal_plane_fin(spacing, D, diffusionTime, 20)(qPoints)

plt.plot(qPoints, trueSignalPoints, color="r")#colors[t], marker=".")
plt.errorbar(qPoints, simulatedSignal, stdsSample, fmt=".", color="g")
plt.legend(["Free diffusion signal", "DIffusion in a sphere","Simulated signal"])
plt.xlabel("q=gamma*G*delta [um-1]")
plt.ylabel("Signal attenuation")
plt.title("graphTitle")
plt.yscale("log")
plt.grid()
plt.show()
"""

diffusionTime = 21
spacing = 16
qPoints = np.linspace(0, 1, 101)[1:]

radius = spacing/2
cylPoints  = Util.getSignal_cylinder_fin(radius, D, diffusionTime, 20, 20, 5)(qPoints)

plt.plot(qPoints, cylPoints)
plt.yscale("log")
plt.show()
