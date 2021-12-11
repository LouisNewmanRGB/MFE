import time
import numpy as np

from SimulationPython.Simulation.Validation import Validation
from SimulationPython.Simulation.Util import Util

if __name__ == '__main__':
    t0 = time.time()
    nRuns = 1
    diffusionTimes = [10] #[2, 5, 25] #ms
    nParts = [38461] #, 100, 1000, 5000, 10000, 100000]
    nStep = 26
    D = 2 #um2/ms
    radius = 8 #um
    T2 = np.inf #no T2 relaxation
    plotGraph = True
    qPoints = np.linspace(0, 0.4, 101)[1:]
    saveFileName = "sphere_conv_29nov2"

    #results = Validation.runParallel(nRuns, Validation.runValidation, [diffusionTimes, (nParts, nStep), "convergence", "sphere_uniform", D, T2, radius], saveFileName=saveFileName)
    #results = Validation.runValidation(nRuns, diffusionTimes, (nParts, nStep), "convergence", "sphere_uniform", D, T2, radius)
    #np.save(Util.getFilePath(saveFileName), results)
    results = np.load(Util.getFilePath(saveFileName)+".npy", allow_pickle=True)
    print(results.shape)

    print("Total time", time.time() - t0)
    plotTitle = "Random walk simulation of diffusion in an impermeable sphere for different diffusion times and numbers of particles\n"\
                "(Number of steps = {n}, {nRuns} run average)".format(n=nStep, nRuns=nRuns)
    Validation.signalPlot(results, plotGraph, plotTitle, qPoints, nRuns, diffusionTimes, (nParts, nStep), "convergence", "sphere_uniform", D, radius, maxStd=3e-3)
    #Validation.averageSignalPlot(results, plotGraph, plotTitle, qPoints, nRuns, diffusionTimes, (nParts, nStep), "convergence", "sphere_uniform", D, radius)
