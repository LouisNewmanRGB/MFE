import time

from Validation import Validation

if __name__ == '__main__':
    t0 = time.time()
    nRuns = 8
    diffusionTimes = [2, 3]#, 10] #ms
    totalSteps = int(1e5)
    nStep = [2, 4, 8, 16]#, 32, 50, 100, 250, 500, 1000] #dividers of 100 000
    #nStep = [6, 7, 8, 9, 10, 12, 16, 20, 35, 50]
    #nStep = [8, 4, 5, 6, 7, 8]
    D = 2 #um2/ms
    radius = 8 #um
    T2 = 1 #irrelevant for this test
    plotHist = False

    errors = Validation.runParallel(nRuns, plotHist, Validation.runSphereComp, [diffusionTimes, nStep, totalSteps, radius, D, T2])
    #errors = Validation.runSphereComp(nRuns, plotHist, diffusionTimes, nStep, totalSteps, radius, D, T2)
    print("Total time", time.time() - t0)
    plotTitle = "Random walk simulation of diffusion in an impermeable sphere for different diffusion times and numbers of steps\n"\
                "particles x steps = {totalSteps}, {nRuns} run average".format(totalSteps=totalSteps, nRuns=nRuns)
    Validation.comparisonFinalPlot(errors, diffusionTimes, nStep, plotTitle)
