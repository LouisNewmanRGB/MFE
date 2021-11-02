import time
import numpy as np

from Validation import Validation
from Util import Util

if __name__ == '__main__':
    t0 = time.time()
    nRuns = 24
    diffusionTimes = [50, 100] #[2, 5, 25] #ms
    nParts = [100, 1000, 5000, 10000]#, 100000]
    nStep = 8
    D = 2 #um2/ms
    radius = 8 #um
    T2 = np.inf #no T2 relaxation
    plotGraph = False
    qPoints = np.linspace(0, 0.4, 101)[1:]
    saveFileName = "sphere_conv_longDT"

    #results = Validation.runParallel(nRuns, Validation.runValidation, [diffusionTimes, (nParts, nStep), "convergence", "sphere_uniform", D, T2, radius], saveFileName=saveFileName)
    #results = Validation.runValidation(nRuns, diffusionTimes, (nParts, nStep), "convergence", "sphere_uniform", D, T2, radius)
    results = np.load(Util.getFilePath(saveFileName)+".npy", allow_pickle=True)
    print(results.shape)

    print("Total time", time.time() - t0)
    plotTitle = "Random walk simulation of diffusion in an impermeable sphere for different diffusion times and numbers of particles\n"\
                "(Number of steps = {n}, {nRuns} run average)".format(n=nStep, nRuns=nRuns)
    #Validation.signalPlot(results, plotGraph, plotTitle, qPoints, nRuns, diffusionTimes, (nParts, nStep), "convergence", "sphere_uniform", D, radius)
    Validation.averageSignalPlot(results, plotGraph, plotTitle, qPoints, nRuns, diffusionTimes, (nParts, nStep), "convergence", "sphere_uniform", D, radius)
