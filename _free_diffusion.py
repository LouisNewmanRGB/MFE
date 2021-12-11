import time

from Validation import Validation

if __name__ == '__main__':
    t0 = time.time()
    nRuns = 25
    diffusionTimes = [1, 20, 100] #ms
    totalSteps = int(1e5)
    nSteps = [2, 4, 8, 16, 32, 50, 100, 250, 500, 1000] #dividers of 100 000
    #nStep = [5, 6, 7, 8, 9, 10, 11, 12, 16, 20]
    D = 2 #um2/ms
    T2 = 1 #irrelevant for this test
    plotHist = False

    results = Validation.runValidation(nRuns, diffusionTimes, (nSteps, totalSteps), "comparison", "free", D, T2)
    print("Total time", time.time() - t0)
    plotTitle = "Random walk simulation of free diffusion for different diffusion times and numbers of steps\n"\
                "particles x steps = {totalSteps}, {nRuns} run average".format(totalSteps=totalSteps, nRuns=nRuns)
    Validation.distributionPlot(results, plotTitle, plotHist, nRuns, diffusionTimes, (nSteps, totalSteps), "comparison", "free", D)
