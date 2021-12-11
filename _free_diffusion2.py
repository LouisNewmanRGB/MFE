import time

from Validation import Validation

if __name__ == '__main__':
    t0 = time.time()
    nRuns = 8
    diffusionTimes = [1, 20, 100] #ms
    nParts = [10, 100, 1000, 10000, 100000]
    nStep = 40
    D = 2 #um2/ms
    T2 = 1 #irrelevant for this test
    plotHist = False

    results = Validation.runValidation(nRuns, diffusionTimes, (nParts, nStep), "convergence", "free", D, T2)
    print("Total time", time.time() - t0)
    plotTitle = "Random walk simulation of free diffusion for different diffusion times and numbers of particles\n"\
                "(Number of steps = {n}, {nRuns} run average)".format(n=nStep, nRuns=nRuns)
    Validation.distributionPlot(results, plotTitle, plotHist, nRuns, diffusionTimes, (nParts, nStep), "convergence", "free", D)
