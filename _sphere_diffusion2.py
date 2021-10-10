import time

from Validation import Validation

if __name__ == '__main__':
    t0 = time.time()
    nRuns = 8
    diffusionTimes = [2, 3, 10] #ms
    nParts = [100, 1000, 10000]#, 100000]
    nStep = 8
    D = 2 #um2/ms
    radius = 8 #um
    T2 = 1 #irrelevant for this test
    plotHist = False

    errors = Validation.runParallel(nRuns, plotHist, Validation.runSphereConv, [diffusionTimes, nStep, nParts, radius, D, T2])
    #errors = Validation.runSphereConv(nRuns, plotHist, diffusionTimes, nStep, nParts, radius, D, T2)
    print("Total time", time.time() - t0)
    plotTitle = "Random walk simulation of diffusion in an impermeable sphere for different diffusion times and numbers of particles\n"\
                "(Number of steps = {n}, {nRuns} run average)".format(n=nStep, nRuns=nRuns)
    Validation.convergenceFinalPlot(errors, diffusionTimes, nParts, plotTitle)
