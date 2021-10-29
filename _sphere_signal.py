import time
import numpy as np

from Validation import Validation
from Util import Util

if __name__ == '__main__':
    t0 = time.time()
    nRuns = 24
    totalSteps = int(1e4)
    diffusionTimes = [3, 10, 50] #ms
    #nStep = [2, 4, 8, 16, 32, 50, 100, 250, 500, 1000] #dividers of 100 000
    nStep = [6, 8, 10, 12, 14, 16, 18, 20, 24, 32]
    D = 2 #um2/ms
    radius = 8 #um
    T2 = np.inf #no T2 relaxation
    plotGraph = False
    qPoints = np.linspace(0, 0.4, 101)[1:]

    results = Validation.runParallel(nRuns, Validation.runSphereSignalComp, [plotGraph, diffusionTimes, nStep, totalSteps, radius, D, T2, qPoints])
    #results = Validation.runSphereSignalComp(nRuns, plotGraph, diffusionTimes, nStep, totalSteps, radius, D, T2, qPoints)
    #results = np.load(Util.getFilePath(Validation.runSphereSignalComp)+".npy")

    errors = results[:, :, 0, :]
    meanDisps = results[:, :, 1, :]
    RMSDs = results[:, :, 2, :]

    print("Total time", time.time() - t0)
    plotTitle = "Random walk simulation of diffusion in an impermeable sphere of radius {radius}um for different diffusion times and numbers of steps\n"\
                "particles x steps = {totalSteps}, {nRuns} run average".format(radius=radius, totalSteps=totalSteps, nRuns=nRuns)
    Validation.signalCompPlot(errors, diffusionTimes, nStep, plotTitle)
    Validation.meanTestCompMD(meanDisps, diffusionTimes, nStep, plotTitle)
    Validation.meanTestCompRMSD(RMSDs, diffusionTimes, plotTitle, nStep, radius, D)
