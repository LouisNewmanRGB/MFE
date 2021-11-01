import time
import numpy as np

from Validation import Validation
from Util import Util

if __name__ == '__main__':
    t0 = time.time()
    nRuns = 8
    totalSteps = int(1e4)
    diffusionTimes = [3, 10, 50] #ms
    nSteps = [2, 4, 8, 16, 32, 50, 100, 250, 500, 1000] #dividers of 100 000
    D = 2 #um2/ms
    nStep = 8
    spacing = 8 #um
    T2 = np.inf #no T2 relaxation
    plotGraph = False
    qPoints = np.linspace(0, 0.4, 101)[1:]
    saveFileName = "planes_comp"

    #results = Validation.runParallel(nRuns, Validation.runValidation, [diffusionTimes, (nSteps, totalSteps), "comparison", "planes", D, T2, spacing], saveFileName)
    results = Validation.runValidation(nRuns, diffusionTimes, (nSteps, totalSteps), "comparison", "planes", D, T2, spacing)
    #results = np.load(Util.getFilePath(saveFileName)+".npy", allow_pickle=True)

    print("Total time", time.time() - t0)
    plotTitle = "MRI signal attenuation between two impermeable planes separated by {spacing}um for different diffusion times\n" \
                "particles x steps = {totalSteps}, {nRuns} run average".format(spacing=spacing, totalSteps=totalSteps, nRuns=nRuns)
    Validation.signalPlot(results, plotGraph, plotTitle, qPoints, nRuns, diffusionTimes, (nSteps, totalSteps), "comparison", "planes", D, spacing)
    Validation.averageSignalPlot(results, plotGraph, plotTitle, qPoints, nRuns, diffusionTimes, (nSteps, totalSteps), "comparison", "planes", D, spacing)
