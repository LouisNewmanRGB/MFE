import time

from SimulationPython.Simulation.Validation import Validation

if __name__ == '__main__':
    t0 = time.time()
    nRuns = 8
    diffusionTimes = [1, 20, 100] #ms
    nParts = [10, 100, 1000, 10000] #, 100000]
    nStep = 8
    D = 2 #um2/ms
    T2 = 1 #irrelevant for this test
    plotHist = False

    errors = Validation.runParallel(nRuns, plotHist, Validation.runFreeConv, [diffusionTimes, nStep, nParts, D, T2])
    #errors = Validation.runFreeConv(nRuns, plotHist, diffusionTimes, nStep, nParts, D, T2)
    print("Total time", time.time() - t0)
    plotTitle = "Random walk simulation of free diffusion for different diffusion times and numbers of particles\n"\
                "(Number of steps = {n}, {nRuns} run average)".format(n=nStep, nRuns=nRuns)
    Validation.convergenceFinalPlot(errors, diffusionTimes, nParts, plotTitle)
