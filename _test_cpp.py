import numpy as np
import matplotlib.pyplot as plt
import scipy.stats

from SimulationResults2 import SimulationResults2
from Util import Util
D = 2
diffusionTime = 50
radius = 8

#data = np.load("./RandomWalkSimulator/results/npyTest.npy")
data = np.load("./RandomWalkSimulator/results/npyTestSphere.npy")
print(data.shape)
print(data[0:7, :])

results = SimulationResults2(data)

"""
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(results.startingPositions[:, 0], results.startingPositions[:, 1], results.startingPositions[:, 2])
#ax.scatter(results.getEndPositions()[:, 0], results.getEndPositions()[:, 1], results.getEndPositions()[:, 2])
#ax.scatter(results.endPositionsVoxel[:, 0], results.endPositionsVoxel[:, 1], results.endPositionsVoxel[:, 2])
plt.show()
"""

qPoints = np.linspace(0, 0.4, 101)[1:]
qVectors = np.array([[q, 0, 0] for q in qPoints])
simulatedSignal, stdsSample, stdsPopulation = results.getSGPSignal(qVectors, includeStd=True)
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


"""
truePDF = Util.getPDF(D, diffusionTime)
trueCDF = Util.getCDF(D, diffusionTime)
pdfPointsX = np.linspace(0, 4 * (6 * D * diffusionTime)**0.5, 500)
pdfPointsY = truePDF(pdfPointsX)
histogramTitle = "TODO"

bw = 2*scipy.stats.iqr(results.getDistances(), rng=(25, 75))/(len(results.getDistances()))**(1/3)
nBins = int((np.max(results.getDistances()) - np.min(results.getDistances()))/bw)
plt.hist(results.getDistances(), bins=nBins, density = True, stacked=True)
plt.plot(pdfPointsX, pdfPointsY, color = 'red')
plt.legend(["Expected probability density function", "Random walk results histogram"])
plt.xlabel("Distance travelled [um]")
plt.ylabel("[um-1]")
plt.title(histogramTitle)
plt.show()
"""
"""
plt.hist(results.getDistances(), bins=50)
print(np.average([d**2 for d in results.getDistances()]))
plt.show()
"""
