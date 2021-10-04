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
    n = nStep[i]
    diffusionTime = diffusionTimes[t]
    timeStep = diffusionTime / n
    nPart = int(totalSteps / n)
    l = (6*D*timeStep)**0.5
    envSize = 5*radius
    env = Environment(T2, D, envSize, envSize, envSize)
    part = [Particle3D(0, 0, 0) for i in range(nPart)]
    sim = Simulation(n, timeStep, part, env, [Sphere(0, 0, 0, T2, D, 0, radius)])
    #sim = Simulation(n, timeStep, part, env)
    sim.run(seed=None, calcData=False)

    #kolmogorov-smirnov
    print("{diffusionTime}ms diffusion time ({t}/{totDt}), {nPart} particles and {n} steps ({i}/{totStep}), run {r}/{nRuns}:"
    .format(diffusionTime=diffusionTime, t=t+1, totDt=len(diffusionTimes), nPart=nPart, n=n, i=i+1, totStep=len(nStep), r=r+1, nRuns=nRuns))
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
        plt.plot(points, distribPoints, color = 'red')
        plt.legend(["Expected probability density function", "Random walk results histogram"])
        plt.xlabel("Distance travelled [um]")
        plt.ylabel("[um-1]")
        plt.title("Probability distribution of the distance travelled by a particle (sphere radius = {radius}um)\n"
                  "Diffusion time = {diffusionTime}ms, Number of particles = {nPart}, Number of steps = {n} ({totalSteps} particles x steps), run {r}/{nRuns}"
                  .format(radius=radius, diffusionTime=diffusionTime, nPart=nPart, n=n, totalSteps=totalSteps, r=r+1, nRuns=nRuns))
        plt.show()

    #numbers
    #print("numer:", np.average(np.power(sim.getDistances(), 2)))
    #print("theor:", 6*D*diffusionTime)

nRuns = 1
diffusionTimes = [10, 3, 10] #ms
totalSteps = int(1e5)
#nStep = [2, 4, 8, 16, 32, 50, 100, 250, 500, 1000] #dividers of 100 000
#nStep = [6, 7, 8, 9, 10, 12, 16, 20, 35, 50]
nStep = [8, 4, 5, 6, 7, 8]
D = 2 #um2/ms
radius = 8 #um
T2 = 1 #irrelevant for this test

errors = np.zeros((len(diffusionTimes), len(nStep), nRuns))
pvalues = np.zeros((len(diffusionTimes), len(nStep), nRuns))

for t in range(len(diffusionTimes)):
    diffusionTime = diffusionTimes[t]
    truePDF = Util.getPDF_sphere(D, diffusionTime, radius, 500, 5)
    trueCDF = Util.getCDF_sphere(D, diffusionTime, radius, 500, 5)

    points = np.linspace(0, 1.1*radius, 500)
    distribPoints = [truePDF(p) for p in points]
    for i in range(len(nStep)):
        for r in range(nRuns):
            runSim(i, t, r, True)

#final plot
print("ERRORS:", errors)
print("PVALUES:", pvalues)
averageErrors = np.average(errors, axis=2)
print("AVERAGE ERRORS:", averageErrors)
colors = cm.rainbow(np.linspace(0, 1, len(diffusionTimes)))
for t in range(len(diffusionTimes)):
    plt.scatter(nStep, averageErrors[t,:], color = colors[t])
plt.xscale("log")
plt.xlabel("Number of steps")
plt.ylabel("Supremum of distances between pseudo-exact and empirical cumulative distribution functions")
plt.legend(["Diffusion time = {diffusionTime}ms".format(diffusionTime=dt) for dt in diffusionTimes])
plt.title("Random walk simulation of diffusion in an impermeable sphere for different diffusion times and numbers of steps\n"
          "particles x steps = {totalSteps}, {nRuns} run average".format(totalSteps=totalSteps, nRuns=nRuns))
plt.show()
