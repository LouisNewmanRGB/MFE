import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy.stats
import time

from Particle3D import Particle3D
from Environment import Environment
from Simulation import Simulation
from Sphere import Sphere
from Util import Util

def runSim(i, t, r, plotHist):
    startTime = time.time()
    nPart = nParts[i]
    diffusionTime = diffusionTimes[t]
    timeStep = diffusionTime / nStep
    l = (6*D*timeStep)**0.5
    envSize = 5*radius
    env = Environment(T2, D, envSize, envSize, envSize)
    #part = [Particle3D(*Util.getRandomDirection()*Util.getRandomQuadratic(radius)) for i in range(nPart)]
    part = [Particle3D(0, 0, 0) for i in range(nPart)]
    for p in part:
        rad = np.linalg.norm(p.getPos())
        if rad > 8:
            print("PROB")
    sim = Simulation(nStep, timeStep, part, env, [Sphere(0, 0, 0, T2, D, 0, radius)])
    sim.run(seed=None, calcData=False)

    #kolmogorov-smirnov
    print("{diffusionTime}ms diffusion time ({t}/{totDt}) and {nPart} particles ({i}/{totStep}), run {r}/{nRuns}:"
    .format(diffusionTime=diffusionTime, t=t+1, totDt=len(diffusionTimes), nPart=nPart, i=i+1, totStep=len(nParts), r=r+1, nRuns=nRuns))
    test = scipy.stats.kstest(sim.getDistances(), trueCDF)
    print("Computation time: {compTime}s".format(compTime = time.time() - startTime))
    print("Kolmogorov-Smirnov test Statistic:", test.statistic)
    print("Kolmogorov-Smirnov test pvalue:", test.pvalue, "\n")
    errors[t, i, r] = test.statistic
    pvalues[t, i, r] = test.pvalue

    #histograms
    if plotHist:
        bw = 2*scipy.stats.iqr(sim.getDistances(), rng=(25, 75))/(len(sim.getDistances()))**(1/3)
        nBins = int((np.max(sim.getDistances()) - np.min(sim.getDistances()))/bw)

        plt.hist(sim.getDistances(), bins=nBins, density = True, stacked=True)
        #plt.hist([np.linalg.norm(pos) for pos in sim.getEndPositions()], bins=nBins, density = True, stacked=True)
        plt.plot(points, distribPoints, color = 'red')
        plt.legend(["Expected probability density function", "Random walk results histogram"])
        plt.xlabel("Distance travelled [um]")
        plt.ylabel("[um-1]")
        plt.title("Probability distribution of the distance travelled by a particle (sphere radius = {radius}um)\n"
                  "Diffusion time = {diffusionTime}ms, Number of particles = {nPart}, Number of steps = {n}, run {r}/{nRuns}"
                  .format(diffusionTime=diffusionTime, nPart=nPart, n=nStep, radius=radius, r=r+1, nRuns=nRuns))
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
        plt.scatter(nParts, averageErrors[t,:], color = colors[t])
    plt.plot(finalPoints, finalPoints**(-0.5), color = "black")

    plt.xscale("log")
    plt.xlabel("Number of particles")
    plt.ylabel("Supremum of distances between exact and empirical cumulative distribution functions")
    plt.legend(["Theoretical distances\n(n^-1/2 convergence rate)"] +
               ["Diffusion time = {diffusionTime}ms".format(diffusionTime=dt) for dt in diffusionTimes])
    plt.title("Random walk simulation of diffusion in an impermeable sphere for different diffusion times and numbers of particles\n"
          "(Number of steps = {n}, {nRuns} run average)".format(n=nStep, nRuns=nRuns))
    plt.show()

nRuns = 5
diffusionTimes = [2, 3, 10] #ms
nParts = [100, 1000, 10000, 100000]
nStep = 8
D = 2 #um2/ms
radius = 8 #um
T2 = 1 #irrelevant for this test

errors = np.zeros((len(diffusionTimes), len(nParts), nRuns))
pvalues = np.zeros((len(diffusionTimes), len(nParts), nRuns))

for t in range(len(diffusionTimes)):
    diffusionTime = diffusionTimes[t]
    truePDF = Util.getPDF_sphere(D, diffusionTime, radius, 500, 5)
    trueCDF = Util.getCDF_sphere(D, diffusionTime, radius, 1000, 5)

    points = np.linspace(0, 1.1*radius, 500)
    distribPoints = [truePDF(p) for p in points]
    for i in range(len(nParts)):
        for r in range(nRuns):
            runSim(i, t, r, False)

#final plot
print("ERRORS:", errors)
print("PVALUES:", pvalues)
averageErrors = np.average(errors, axis=2)
print("AVERAGE ERRORS:", averageErrors)
finalPlot(logScale=False)
finalPlot(logScale=True)
