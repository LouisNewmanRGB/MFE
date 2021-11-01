import time
import numpy as np

from Validation import Validation
from Util import Util

if __name__ == '__main__':
    t0 = time.time()
    nRuns = 24
    totalSteps = int(1e4)
    diffusionTimes = [3, 10, 50] #ms
    nSteps = [2, 4, 8, 16, 32, 50, 100, 250, 500, 1000] #dividers of 100 000
    #nSteps = [6, 8, 10, 12, 14, 16, 18, 20, 24, 32]
    D = 2 #um2/ms
    radius = 8 #um
    T2 = np.inf #no T2 relaxation
    plotGraph = False
    qPoints = np.linspace(0, 0.4, 101)[1:]
    saveFileName = "sphere_signal_comp"

    #results = Validation.runParallel(nRuns, Validation.runValidation, [diffusionTimes, (nSteps, totalSteps), "comparison", "sphere_uniform", D, T2, radius], saveFileName=saveFileName)
    #results = Validation.runValidation(nRuns, diffusionTimes, (nSteps, totalSteps), "comparison", "sphere_uniform", D, T2, radius)
    results = np.load(Util.getFilePath(saveFileName)+".npy", allow_pickle=True)

    print("Total time", time.time() - t0)
    plotTitle = "Random walk simulation of diffusion in an impermeable sphere of radius {radius}um for different diffusion times and numbers of steps\n"\
                "particles x steps = {totalSteps}, {nRuns} run average".format(radius=radius, totalSteps=totalSteps, nRuns=nRuns)
    #Validation.signalPlot(results, plotGraph, plotTitle, qPoints, nRuns, diffusionTimes, (nSteps, totalSteps), "comparison", "sphere_uniform", D, radius)
    Validation.averageSignalPlot(results, plotGraph, plotTitle, qPoints, nRuns, diffusionTimes, (nSteps, totalSteps), "comparison", "sphere_uniform", D, radius)
