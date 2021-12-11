import time

from Validation import Validation
#from SimulationPython.Simulation.Validation import Validation

if __name__ == '__main__':
    t0 = time.time()
    nRuns = 10
    diffusionTimes = [2, 3, 50] #ms
    nParts = [100, 1000, 10000, 100000]
    nStep = 8
    D = 2 #um2/ms
    radius = 8 #um
    T2 = 1 #irrelevant for this test
    plotHist = False

    results = Validation.runValidation(nRuns, diffusionTimes, (nParts, nStep), "convergence", "sphere_center", D, T2, parameter=radius)
    print("Total time", time.time() - t0)
    plotTitle = "Random walk simulation of diffusion in an impermeable sphere for different diffusion times and numbers of particles\n"\
                "(Number of steps = {n}, {nRuns} run average)".format(n=nStep, nRuns=nRuns)
    Validation.distributionPlot(results, plotTitle, plotHist, nRuns, diffusionTimes, (nParts, nStep), "convergence", "sphere_center", D, parameter=radius)
