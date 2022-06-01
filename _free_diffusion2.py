import time

from Validation import Validation

if __name__ == '__main__':
    t0 = time.time()
    nRuns = 50
    diffusionTimes = [1, 20, 100] #ms
    nParts = [10, 100, 1000, 10000, 100000]
    nStep = 8
    D = 2 #um2/ms
    T2 = 1 #irrelevant for this test
    plotHist = False

    results = Validation.runValidation(nRuns, diffusionTimes, (nParts, nStep), "convergence", "free", D, T2)
    print("Total time", time.time() - t0)
    plotTitle = "Monte Carlo Simulation of Free Diffusion\n"\
                "(Number of steps = {n}, {nRuns} run average)".format(n=nStep, nRuns=nRuns)
    Validation.distributionPlot_normal(results, plotTitle, plotHist, nRuns, diffusionTimes, (nParts, nStep), "convergence", "free", D)
