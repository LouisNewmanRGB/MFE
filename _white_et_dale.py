import numpy as np
import matplotlib.pyplot as plt
import time
import scipy.stats

from Environment import Environment
from Sphere import Sphere
from Particle3DPiecewiseGradient import Particle3DPiecewiseGradient
from Simulation import Simulation
from Util import Util
from Validation import Validation
from SimulationResults import SimulationResults
from StejskalTanner import StejskalTanner

#simulation parameters
diffusionTimes = [12, 18, 60] #[12, 18, 24, 30, 36, 42, 48, 54, 60] #[12] for white_et_dale_permeability and white_et_dale_permeability2
nPart = int(1e3) #1e4 for white_et_dale_permeability 1e5 for white_et_dale_permeability2
nucleusFractions = [0.8] #[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8] #volume fraction #[0.7] for white_et_dale_permeability2
TE = 140 #ms
nStep = 195 #8
bValue = 4 #ms/um2
saveFileName = "white_et_dale_195TS" #for white_et_dale_195TS: diffusionTimes = [12, 18, 60] and nucleusFractions = [0.8]
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

if load:
    results = np.load(Util.getFilePath(saveFileName)+".npy", allow_pickle=True)
else:
    results = np.empty((len(diffusionTimes), len(nucleusFractions)), dtype=SimulationResults)
    for t in range(len(diffusionTimes)):
        Delta = diffusionTimes[t]
        timeStep = TE / nStep
        delta = Delta
        delay = TE/2 - Delta
        seq = StejskalTanner(delay, TE, delta, Delta)
        for f in range(len(nucleusFractions)):
            startTime = time.time()
            fraction = nucleusFractions[f]
            nucleusRadius = fraction**(1/3) * cellRadius
            nucleus = Sphere(0, 0, 0, T2N, diffusivityN, permeabilityNC, nucleusRadius)
            particles = [Particle3DPiecewiseGradient(*Util.getRandomDirection()*Util.getRandomQuadratic(cellRadius), seq) for i in range(nPart)]
            sim = Simulation(nStep, timeStep, particles, env, [nucleus, cytoplasm])
            sim.run(calcData=False)
            #sim.plot(positionType=False)
            results[t, f] = sim.getResults()
            print( "{diffusionTime}ms diffusion time ({t}/{totDt}), {frac} nucleus fraction ({f}/{nFrac}):"\
                        .format(diffusionTime=Delta, t=t+1, totDt=len(diffusionTimes), frac = fraction, f=f+1, nFrac = len(nucleusFractions)) )
            print("computation time:", time.time() - startTime, "s")
    np.save(Util.getFilePath(saveFileName), results)

nStepMin1 = Util.getMinNStep(TE, permeabilityNC, diffusivityN, diffusivityC)
print("Minimum number of time steps required:", nStepMin1)

f = 0
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
        signals[f], stdsSample[f], stdsPopulation[f] = (results[t, f]).getFiniteGradientSignalBValue(bValue, includeStd=True)
    plt.errorbar(nucleusFractions, signals, stdsSample, fmt=".")
plt.title("Simulation of MRI signal in cells with different nuclear volume fractions\n"
          "T2cytosol = {T2C}ms, T2nucleus = {T2N}ms".format(T2C=T2C, T2N=T2N))
plt.ylabel("MRI signal attenuation")
plt.xlabel("Nuclear volume fraction")
plt.legend(["Diffusion time = {diffusionTime}ms".format(diffusionTime=dt) for dt in diffusionTimes])
plt.grid()
plt.show()


