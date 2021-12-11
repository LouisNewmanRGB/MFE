import time
import numpy as np

from SimulationPython.Simulation.Validation import Validation

if __name__ == '__main__':
    t0 = time.time()
    nRuns = 24
    diffusionTimes = [3, 10, 50] #ms
    nParts = [100, 1000, 10000]#, 100000]
    nStep = 8

    D = 2 #um2/ms
    spacing = 8 #um
    T2 = np.inf #no T2 relaxation
    plotGraph = False
    qPoints = np.linspace(0, 0.4, 101)[1:]
    saveFileName = "planes_conv"

    results = Validation.runParallel(nRuns, Validation.runValidation, [diffusionTimes, (nParts, nStep), "convergence", "planes", D, T2, spacing], saveFileName=saveFileName)
    #results = Validation.runValidation(nRuns, diffusionTimes, (nParts, nStep), "convergence", "planes", D, T2, spacing)
    #results = np.load(Util.getFilePath(saveFileName)+".npy", allow_pickle=True)
    print(results.shape)

    print("Total time", time.time() - t0)
    plotTitle = "Random walk simulation of diffusion between two impermeable planes for different diffusion times and numbers of particles\n"\
                "(Number of steps = {n}, {nRuns} run average)".format(n=nStep, nRuns=nRuns)
    #Validation.signalPlot(results, plotGraph, plotTitle, qPoints, nRuns, diffusionTimes, (nParts, nStep), "convergence", "planes", D, spacing)
    Validation.averageSignalPlot(results, plotGraph, plotTitle, qPoints, nRuns, diffusionTimes, (nParts, nStep), "convergence", "planes", D, spacing)
