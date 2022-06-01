import numpy as np
import matplotlib.pyplot as plt

from SimulationCppResults import SimulationCppResults
from SimulationCpp import Simulation
from Util import Util

#diffusionTimes = list(np.arange(1, 50, 2)) + list(np.arange(52, 76, 4))
diffusionTimes = [10]
nStep = 26
nPart = int(308461)
D = 2
T2 = np.inf
radius = 8
load = False
plotIntermediary = True

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

    qPoints = np.linspace(0, 0.8, 101)[1:]
    qVectors = np.array([[q, 0, 0] for q in qPoints])

    RMSEs = np.empty((len(diffusionTimes)))
    for t in range(len(diffusionTimes)):
        diffusionTime = diffusionTimes[t]
        trueSignalPoints = Util.getSignal_sphere_fin(radius, D, diffusionTime, 20, 20, 10)(qPoints)
        simResults = results[t]
        simulatedSignal, stdsSample, stdsPopulation = simResults.getSGPSignal(qVectors, includeStd=True)
        RMSEs[t] = ( np.average((trueSignalPoints - simulatedSignal)**2) )**0.5

        if plotIntermediary:
            plt.plot(qPoints, trueSignalPoints, color="r")#colors[t], marker=".")
            plt.errorbar(qPoints, simulatedSignal, stdsSample, fmt=".", color="g")
            plt.legend(["Theoretical signal", "Simulated signal", "Simulated signal (real part)"])
            plt.xlabel("q=gamma*G*delta [um-1]")
            plt.ylabel("Signal attenuation")
            plt.title("Monte Carlo Simulation of the SGP Signal for Diffusion in an Impermeable Sphere")
            plt.yscale("log")
            plt.grid()
            plt.show()
    #np.save("numpy saves/sphere_signalcpp24.npy", RMSEs)
else:
    RMSE8 = np.load("numpy saves/sphere_signalcpp8.npy", allow_pickle=True)
    RMSE16 = np.load("numpy saves/sphere_signalcpp16.npy", allow_pickle=True)
    RMSE24 = np.load("numpy saves/sphere_signalcpp24.npy", allow_pickle=True)

plt.scatter(diffusionTimes, RMSE8)
plt.scatter(diffusionTimes, RMSE16)
#plt.scatter(diffusionTimes, RMSE24)
plt.xlabel("Diffusion time [ms]")
plt.ylabel("RMSE")
plt.legend(["8 steps", "16 steps"])
plt.title("Theoretical SGP Signal Attenuation for Diffusion in an Impermeable Sphere")
plt.show()

