import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy.stats
import time

from Particle3D import Particle3D
from Environment import Environment
from Simulation import Simulation
from Util import Util

def runSim(i, t, plotHist):
    startTime = time.time()
    nPart = nParts[i]
    diffusionTime = diffusionTimes[t]
    timeStep = diffusionTime / nStep
    l = (6*D*timeStep)**0.5
    envSize = 10*l
    env = Environment(T2, D, envSize, envSize, envSize)
    part = [Particle3D(Util.getRandomU(envSize),Util.getRandomU(envSize),Util.getRandomU(envSize)) for i in range(nPart)]
    sim = Simulation(nStep, timeStep, part, env)
    sim.run(seed=None, calcData=False)

    #kolmogorov-smirnov
    test = scipy.stats.kstest(sim.getDistances(), trueCDF)
    print("{diffusionTime}ms diffusion time ({t}/{totDt}) and {nPart} particles ({i}/{totStep}):"
    .format(diffusionTime=diffusionTime, t=t+1, totDt=len(diffusionTimes), nPart=nPart, i=i+1, totStep=len(nParts)))
    print("Computation time: {compTime}s".format(compTime = time.time() - startTime))
    print("Kolmogorov-Smirnov test Statistic:", test.statistic)
    print("Kolmogorov-Smirnov test pvalue:", test.pvalue, "\n")
    errors[t, i] = test.statistic
    pvalues[t, i] = test.pvalue

    #histograms
    if plotHist:
        bw = 2*scipy.stats.iqr(sim.getDistances(), rng=(25, 75))/(len(sim.getDistances()))**(1/3)
        nBins = int((np.max(sim.getDistances()) - np.min(sim.getDistances()))/bw)
        plt.hist(sim.getDistances(), bins=nBins, density = True, stacked=True)
        plt.plot(points, distribPoints, color = 'red')
        plt.legend(["Expected probability density function", "Random walk results histogram"])
        plt.xlabel("Distance travelled [mm]")
        plt.ylabel("[mm-1]")
        plt.title("Probability distribution of the distance travelled by a particle\n"
                  "Diffusion time = {diffusionTime}ms, Number of particles = {nPart}, Number of steps = {n}"
                  .format(diffusionTime=diffusionTime, nPart=nPart, n=nStep,))
        plt.show()

    #numbers
    #print("numer:", np.average(np.power(sim.getDistances(), 2)))
    #print("theor:", 6*D*diffusionTime)

def finalPlot(logScale=False):
    if logScale:
        plt.yscale("log")
        finalPoints = np.array([nParts[0]/10, nParts[-1]*10])
    else:
        finalPoints = 10**np.linspace(np.log10(nParts[0]/1.5), np.log10(nParts[-1]*5), 500)

    colors = cm.rainbow(np.linspace(0, 1, len(diffusionTimes)))
    for t in range(len(diffusionTimes)):
        plt.scatter(nParts, errors[t,:], color = colors[t])
    plt.plot(finalPoints, finalPoints**(-0.5), color = "black")

    plt.xscale("log")
    plt.xlabel("Number of particles")
    plt.ylabel("Supremum of distances between exact and empirical cumulative distribution functions")
    plt.legend(["Theoretical distances\n(n^-1/2 convergence rate)"] +
               ["Diffusion time = {diffusionTime}ms".format(diffusionTime=dt) for dt in diffusionTimes])
    plt.title("Random walk simulation of free diffusion for different diffusion times and numbers of particles\n"
          "(Number of steps = {n})".format(n=nStep))
    plt.show()

diffusionTimes = [1, 20, 100] #ms
nParts = [10, 100, 1000, 10000]#, 100000]
nStep = 8
D = 2e-3/1000 #mm2/ms
T2 = 1 #irrelevant for this test

errors = np.zeros((len(diffusionTimes), len(nParts)))
pvalues = np.zeros((len(diffusionTimes), len(nParts)))

for t in range(len(diffusionTimes)):
    diffusionTime = diffusionTimes[t]
    truePDF = Util.getPDF(D, diffusionTime)
    trueCDF = Util.getCDF(D, diffusionTime)

    points = np.linspace(0, 4 * (6 * D * diffusionTime)**0.5, 500)
    distribPoints = truePDF(points)
    for i in range(len(nParts)):
        tim = time.time()
        runSim(i, t, False)
        #print("time:", time.time() - tim)

#final plot
finalPlot(logScale=False)
finalPlot(logScale=True)
