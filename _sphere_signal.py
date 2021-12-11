import numpy as np
import matplotlib.pyplot as plt

from SimulationCppResults import SimulationCppResults
from SimulationCpp import Simulation
from Util import Util

diffusionTimes = np.array([1, 2, 3] + list(np.linspace(1, 25, 25)*4) )
print(diffusionTimes)
nStep = 8
nPart = int(1e6)
D = 2
T2 = np.inf
radius = 8
load = True

if not(load):
    results = np.empty((len(diffusionTimes)), dtype=SimulationCppResults)
    for t in range(len(diffusionTimes)):
        diffusionTime = diffusionTimes[t]
        timeStep = diffusionTime/nStep
        envSize = radius*5
        sim = Simulation(nStep, timeStep, T2, D, envSize, envSize, envSize)
        sim.addSphere(0, 0, 0, T2, D, 0, radius)
        sim.createStartingPositions(nPart, True)
        sim.createSequenceSGP()
        sim.run()
        results[t] = SimulationCppResults(sim.getResults(), nPart)
        print("Finished sim", t+1)
    np.save("numpy saves/sphere_signalcpp.npy", results)
else:
    np.load("numpy saves/sphere_signalcpp.npy")
qPoints = np.linspace(0, 0.4, 101)[1:]
qVectors = np.array([[q, 0, 0] for q in qPoints])
plotIntermediary = False



RMSEs = np.empty((len(diffusionTimes)))
STDs = np.empty((len(diffusionTimes)))
for t in range(len(diffusionTimes)):
    diffusionTime = diffusionTimes[t]
    trueSignalPoints = Util.getSignal_sphere_fin(radius, D, diffusionTime, 20, 20, 10)(qPoints)
    simResults = results[t]
    simulatedSignal, stdsSample, stdsPopulation = simResults.getSGPSignal(qVectors, includeStd=True)
    RMSEs[t] = ( np.average((trueSignalPoints - simulatedSignal)**2) )**0.5
    #STDs[t] = ( np.std((trueSignalPoints - simulatedSignal)**2) )**0.5 #???????????????????????????????????????????/
    STDs[t] = (np.average(stdsSample**2))**0.5 #???????????????????????????????????????????

    if plotIntermediary:
        plt.plot(qPoints, trueSignalPoints, color="r")#colors[t], marker=".")
        plt.errorbar(qPoints, simulatedSignal, stdsSample, fmt=".", color="g")
        plt.legend(["Theoretical signal", "Simulated signal", "Simulated signal (real part)"])
        plt.xlabel("q=gamma*G*delta [um-1]")
        plt.ylabel("SGP Signal attenuation")
        plt.title("graphTitle")
        plt.yscale("log")
        plt.grid()
        plt.show()

plt.errorbar(diffusionTimes, RMSEs, STDs, fmt="o")
plt.ylabel("Diffusion time [ms]")
plt.xlabel("RMSE between theoretical and simulated SGP signals")
plt.show()

