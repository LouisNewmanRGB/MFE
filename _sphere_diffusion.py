import time

from Validation import Validation

if __name__ == '__main__':
    t0 = time.time()
    nRuns = 50
    diffusionTimes = [2, 5, 10] #ms
    totalSteps = int(1e5)
    nSteps = [2, 4, 8, 16, 32, 50, 100, 250, 500, 1000] #dividers of 100 000
    D = 2 #um2/ms
    radius = 8 #um
    T2 = 1 #irrelevant for this test
    plotHist = False

    results = Validation.runValidation(nRuns, diffusionTimes, (nSteps, totalSteps), "comparison", "sphere_center", D, T2, parameter=radius)
    print("Total time", time.time() - t0)
    plotTitle = "Monte Carlo Simulation of Diffusion in an Impermeable Sphere\n"\
                "particles x steps = {totalSteps}, {nRuns} run average".format(totalSteps=totalSteps, nRuns=nRuns)
    Validation.distributionPlot(results, plotTitle, plotHist, nRuns, diffusionTimes, (nSteps, totalSteps), "comparison", "sphere_center", D, parameter=radius)
