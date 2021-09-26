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
    n = nStep[i]
    diffusionTime = diffusionTimes[t]
    timeStep = diffusionTime / n
    nPart = int(totalSteps / n)
    l = (6*D*timeStep)**0.5
    envSize = 10*l
    env = Environment(T2, D, envSize, envSize, envSize)
    part = [Particle3D(Util.getRandomU(envSize),Util.getRandomU(envSize),Util.getRandomU(envSize)) for i in range(nPart)]
    sim = Simulation(n, timeStep, part, env)
    sim.run(seed=None, calcData=False)

    #kolmogorov-smirnov
    test = scipy.stats.kstest(sim.getDistances(), trueCDF)
    print("{diffusionTime}ms diffusion time ({t}/{totDt}), {nPart} particles and {n} steps ({i}/{totStep}):"
    .format(diffusionTime=diffusionTime, t=t+1, totDt=len(diffusionTimes), nPart=nPart, n=n, i=i+1, totStep=len(nStep)))
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
                  "Diffusion time = {diffusionTime}ms, Number of particles = {nPart}, Number of steps = {n}, ({totalSteps} particles x steps)"
                  .format(diffusionTime=diffusionTime, nPart=nPart, n=n, totalSteps=totalSteps))
        plt.show()

    #numbers
    #print("numer:", np.average(np.power(sim.getDistances(), 2)))
    #print("theor:", 6*D*diffusionTime)

diffusionTimes = [1, 20, 100] #ms
totalSteps = int(1e5)
nStep = [2, 4, 8, 16, 32, 50, 100, 250, 500, 1000] #dividers of 100 000
D = 2e-3/1000 #mm2/ms
T2 = 1 #irrelevant for this test

errors = np.zeros((len(diffusionTimes), len(nStep)))
pvalues = np.zeros((len(diffusionTimes), len(nStep)))

for t in range(len(diffusionTimes)):
    diffusionTime = diffusionTimes[t]
    truePDF = Util.getPDF(D, diffusionTime)
    trueCDF = Util.getCDF(D, diffusionTime)

    points = np.linspace(0, 4 * (6 * D * diffusionTime)**0.5, 500)
    distribPoints = truePDF(points)
    for i in range(len(nStep)):
        tim = time.time()
        runSim(i, t, False)
        #print("time:", time.time() - tim)

#final plot
print("ERRORS:", errors)
print("PVALUES:", pvalues)
colors = cm.rainbow(np.linspace(0, 1, len(diffusionTimes)))
for t in range(len(diffusionTimes)):
    plt.scatter(nStep, errors[t,:], color = colors[t])
plt.xscale("log")
plt.xlabel("Number of steps")
plt.ylabel("Supremum of distances between exact and empirical cumulative distribution functions")
plt.legend(["Diffusion time = {diffusionTime}ms".format(diffusionTime=dt) for dt in diffusionTimes])
plt.title("Random walk simulation of free diffusion for different diffusion times and numbers of steps\n"
          "particles x steps = {totalSteps}".format(totalSteps=totalSteps))
plt.show()
