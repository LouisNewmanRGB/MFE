import time
import numpy as np

from Validation import Validation

if __name__ == '__main__':
    t0 = time.time()
    nRuns = 1
    diffusionTimes = [2, 3, 10] #ms
    #nParts = [10000, 100, 1000, 10000]#, 100000]
    nPart = 10000
    D = 2 #um2/ms
    #magicl2 = 3
    #nStep = int(np.rint(6*D*diffusionTimes[0]/magicl2))
    nStep = 8
    radius = 8 #um
    T2 = np.inf #no T2 relaxation
    plotHist = True
    qPoints = np.linspace(0, 1.5, 101)[1:]

    #errors = Validation.runParallel(nRuns, plotHist, Validation.runSphereConv, [diffusionTimes, nStep, nParts, radius, D, T2])
    signals = Validation.runSphereSignalConv(nRuns, plotHist, diffusionTimes, nStep, nPart, radius, D, T2, qPoints)
    print("Total time", time.time() - t0)
    plotTitle = "MRI signal attenuation in an impermeable sphere of radius {radius}um for different diffusion times\n" \
                "({nPart} particles, {nStep} time steps, {nRuns} run average)".format(radius = radius, nPart=nPart, nStep=nStep, nRuns=nRuns)
    Validation.convSignalFinalPlot(signals, diffusionTimes, radius, D, qPoints, plotTitle)
