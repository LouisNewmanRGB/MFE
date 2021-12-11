import numpy as np

data = np.load("./RandomWalkSimulator/results/white_et_dale_cpp2.npy")
print(data.shape)
print(data[0:7, :])

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats

from SimulationPython.Simulation.Environment import Environment
from SimulationPython.Simulation.Sphere import Sphere
from SimulationPython.Simulation.Util import Util
from SimulationCppResults import SimulationResults2

#simulation parameters
diffusionTimes = [12, 18, 24, 30, 36, 42, 48, 54, 60] #[12] for white_et_dale_permeability and white_et_dale_permeability2
nPart = int(1e3) #1e4 for white_et_dale_permeability 1e5 for white_et_dale_permeability2
nucleusFractions = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8] #volume fraction #[0.7] for white_et_dale_permeability2
TE = 140 #ms
nStep = 195 #8
bValue = 4 #ms/um2
saveFileName = "./RandomWalkSimulator/results/white_et_dale_cpp.npy" #for white_et_dale_195TS: diffusionTimes = [12, 18, 60] and nucleusFractions = [0.8]
load = False

#environment parameters
permeabilityCE = 0
permeabilityNC = 0.5 #0.1 #0.5 #um/ms
diffusivityN = 1.31 #um2/ms
diffusivityC = 0.48 #um2/ms
diffusivityExtra = 1.82 #um2/ms
T2N = 63.29 #ms
T2C = 23.87 #ms
T2Extra = 150 #ms
cellRadius = 5 #um

envSize = 5*cellRadius
env = Environment(T2Extra, diffusivityExtra, envSize, envSize, envSize)
cytoplasm = Sphere(0, 0, 0, T2C, diffusivityC, permeabilityCE, cellRadius)

results = np.empty((len(diffusionTimes), len(nucleusFractions)), dtype=SimulationResults2)
data = np.load(saveFileName)
for t in range(len(diffusionTimes)):
    for f in range(len(nucleusFractions)):
        results[t, f] = SimulationResults2(data[t, f, :, :])

nStepMin1 = Util.getMinNStep(TE, permeabilityNC, diffusivityN, diffusivityC)
print("Minimum number of time steps required:", nStepMin1)

f = 2
fraction = nucleusFractions[f]
pdf = Util.getPDF_sphere_standing(cellRadius)
points = np.linspace(0, cellRadius*1.1, 500)
plt.plot(points, pdf(points), color="r")
interfacePoint = fraction**(1/3) * cellRadius
plt.vlines(interfacePoint, ymin=0, ymax=pdf(interfacePoint), color="r")
data = np.linalg.norm(results[-1, f].getEndPositions(), axis=1)
bw = 2*scipy.stats.iqr(data, rng=(25, 75))/(len(data))**(1/3)
nBins = int((np.max(data) - np.min(data))/bw)
plt.hist(data, bins=nBins, density = True, stacked=True)
plt.grid()
plt.show()

for t in range(len(diffusionTimes)):
    diffusionTime = diffusionTimes[t]
    signals = np.empty(len(nucleusFractions))
    stdsSample = np.empty(len(nucleusFractions))
    stdsPopulation = np.empty(len(nucleusFractions))
    for f in range(len(nucleusFractions)):
        signals[f], stdsSample[f], stdsPopulation[f] = (results[t, f]).getFiniteGradientSignalBValue(bValue, diffusionTime, diffusionTime, includeStd=True)
    plt.errorbar(nucleusFractions, signals, stdsSample, fmt=".")
plt.title("Simulation of MRI signal in cells with different nuclear volume fractions\n"
          "T2cytosol = {T2C}ms, T2nucleus = {T2N}ms".format(T2C=T2C, T2N=T2N))
plt.ylabel("MRI signal attenuation")
plt.xlabel("Nuclear volume fraction")
plt.legend(["Diffusion time = {diffusionTime}ms".format(diffusionTime=dt) for dt in diffusionTimes])
plt.grid()
plt.show()


