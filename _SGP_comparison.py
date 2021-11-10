import numpy as np
import matplotlib.pyplot as plt
import time

from StejskalTanner import StejskalTanner
from SimulationResults import SimulationResults
from Environment import Environment
from Simulation import Simulation
from Sphere import Sphere
from Particle3DPiecewiseGradient import Particle3DPiecewiseGradient
from Util import Util

nPart = 10000
nStep = 8

Delta = 50 #ms
deltaList = [10, 5, 1] #ms
T2 = np.inf
D = 2
radius = 8 #um

qPoints = np.linspace(0, 0.4, 101)[1:]

saveFileName = "SGP_comparison"
load = True
plotSequences = True

if load:
    results = np.load(Util.getFilePath(saveFileName)+".npy", allow_pickle=True)
else:
    results = np.empty(len(deltaList), dtype=SimulationResults)
    for d in range(len(deltaList)):
        delta = deltaList[d]
        TE = 2*Delta
        t = 0
        seq = StejskalTanner(t, TE, delta, Delta)
        if plotSequences:
            seq.plot()
        startTime = time.time()
        timeStep = (Delta + delta)/nStep
        part = [Particle3DPiecewiseGradient(*Util.getRandomDirection()*Util.getRandomQuadratic(radius), seq) for i in range(nPart)]
        envSize = 5*radius
        compartments = [Sphere(0, 0, 0, T2, D, 0, radius)]
        env = Environment(T2, D, envSize, envSize, envSize)
        sim = Simulation(nStep, timeStep, part, env, compartments)
        sim.run()
        results[d] = sim.getResults()
        print("delta = {delta} ({d}/{deltaNumber}):".format(delta=delta, d=d+1, deltaNumber=len(deltaList)))
        print("computation time:", time.time() - startTime, "s", "\n")
    np.save(Util.getFilePath(saveFileName), results)

print("plotting")
theoretical = Util.getSignal_sphere_fin(radius, D, Delta, 50, 50, 10)
theoreticalPoints = theoretical(qPoints)
plt.plot(qPoints, theoreticalPoints, color="r")
for d in range(len(deltaList)):
    delta = deltaList[d]
    gammaGPoints = qPoints/delta
    signal = np.empty(len(qPoints))
    populationStds = np.empty(len(qPoints))
    sampleStds = np.empty(len(qPoints))
    for g in range(len(gammaGPoints)):
        gammaG = gammaGPoints[g]
        signal[g], sampleStds[g], populationStds[g] = results[d].getFiniteGradientSignalGammaG(gammaG, includeStd=True)
    plt.errorbar(qPoints, signal, sampleStds, fmt=".")
plt.title("Simulation of Stejskal-Tanner sequences in a sphere\n"
          "radius={radius}mm, {nPart} particles, {nStep} time steps, {Delta}ms diffusion time".
          format(radius=radius, nPart=nPart, nStep=nStep, Delta=Delta))
plt.xlabel("q = gamma*G*delta [um-1]")
plt.ylabel("MRI signal attenuation")
plt.legend(["Theoretical SGP signal"] + ["delta = {delta}ms".format(delta=delta) for delta in deltaList])
plt.grid()
plt.show()
