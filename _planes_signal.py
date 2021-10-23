import time
import numpy as np

from Validation import Validation
from Util import Util

if __name__ == '__main__':
    t0 = time.time()
    nRuns = 1
    diffusionTimes = [10] #, 3, 10] #ms
    #nParts = [10000, 100, 1000, 10000]#, 100000]
    nPart = 10000
    D = 2 #um2/ms
    #magicl2 = 3
    #nStep = int(np.rint(6*D*diffusionTimes[0]/magicl2))
    nStep = 8
    spacing = 12 #um
    T2 = np.inf #no T2 relaxation
    #plotHistType = "positions_norm"
    #plotHistType = "positions_x"
    plotHistType = "displacements_x"
    qPoints = np.linspace(0, 1.5, 101)[1:]

    #errors = Validation.runParallel(nRuns, plotHist, Validation.runSphereConv, [diffusionTimes, nStep, nParts, radius, D, T2])
    signals = Validation.runPlanesSignalConv(nRuns, plotHistType, diffusionTimes, nStep, nPart, spacing, D, T2, qPoints)
    print("Total time", time.time() - t0)
    plotTitle = "MRI signal attenuation between two impermeable planes separated by {spacing}um for different diffusion times\n" \
                "({nPart} particles, {nStep} time steps, {nRuns} run average)".format(spacing = spacing, nPart=nPart, nStep=nStep, nRuns=nRuns)
    theoreticalSignals = [Util.getSignal_plane_fin(spacing, D, dt, 500) for dt in diffusionTimes] + [Util.getSignal_plane_inf(spacing)]
    #theoreticalSignals = [Util.getSignal_sphere_fin(spacing, D, dt, 10, 10) for dt in diffusionTimes] + [Util.getSignal_sphere_inf(spacing)]
    Validation.convSignalFinalPlot(signals, diffusionTimes, D, qPoints, plotTitle, theoreticalSignals)
