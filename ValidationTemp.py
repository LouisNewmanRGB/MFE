import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import time

from Particle3D import Particle3D
from Environment import Environment
from Simulation import Simulation
from Sphere import Sphere
from Planes import Planes
from Util import Util

class ValidationTemp():

    def run1SimSignal(sim, printMessage, histogramTitle, trueCDF, pdfPointsX, pdfPointsY, plotHistType, qPoints):
        startTime = time.time()
        sim.run(seed=None, calcData=False)
        sim = sim.getResults()
        print("Computation time: {compTime}s".format(compTime = time.time() - startTime))
        print(printMessage)

        if plotHistType != None:
            if plotHistType == "positions_x":
                positions = sim.getEndPositions()[:, 0] #Xs
            elif plotHistType == "positions_norm":
                positions = np.linalg.norm(sim.getEndPositions(), axis=1)
            elif plotHistType == "displacements_x":
                positions = sim.getDisplacements()[:, 0]
            else:
                print("Error: invalid data requested for histogram")
            bw = 2*scipy.stats.iqr(positions, rng=(25, 75))/(len(positions))**(1/3)
            nBins = int((np.max(positions) - np.min(positions))/bw)
            plt.hist(positions, bins=nBins, density = True, stacked=True)
            plt.plot(pdfPointsX, pdfPointsY, color = 'red')
            plt.legend(["Expected probability density function", "Random walk results histogram"])
            plt.xlabel("Distance travelled [um]")
            plt.ylabel("[um-1]")
            plt.title(histogramTitle)
            plt.grid()
            plt.show()

        #signal computations:
        print("Average displacement vector:", sim.getAvgDisplacements())
        qVectors = np.array([[q, 0, 0] for q in qPoints])
        simulatedSignal = sim.getSGPSignal(qVectors, real=False)
        simulatedSignalReal = sim.getSGPSignal(qVectors, real=True)
        return np.array([simulatedSignal, simulatedSignalReal])

    def convSignalFinalPlot(signals, diffusionTimes, D, qPoints, plotTitle, theoreticalSignals):
        averageSignals = np.average(signals, axis=1)
        #colors = cm.rainbow(np.linspace(0, 1, len(diffusionTimes)))
        theoreticalInfinite = theoreticalSignals[-1]
        plt.plot(qPoints, theoreticalInfinite(qPoints), color="black")
        for t in range(len(diffusionTimes)):
            dt = diffusionTimes[t]
            theoreticalSignal = theoreticalSignals[t]
            simulatedSignal = averageSignals[t, 0, :]
            simulatedSignalReal = averageSignals[t, 1, :]
            plt.plot(qPoints, theoreticalSignal(qPoints), color="r")#colors[t], marker=".")
            plt.plot(qPoints, simulatedSignal, color="g")#colors[t], marker="o")
            plt.plot(qPoints, simulatedSignalReal, color="b")#colors[t], marker="v")
        plt.legend(["Theoretical signal for an infinite diffusion time", "Theoretical signal", "Simulated signal", "Simulated signal (2x real part)"])
        plt.xlabel("q=gamma*G*delta [um-1]")
        plt.ylabel("Signal attenuation")
        plt.title(plotTitle)
        plt.yscale("log")
        plt.grid()
        plt.show()

    def runSphereSignalConv(nRuns, plotHistType, diffusionTimes, nStep, nPart, radius, D, T2, qPoints):
        signals = np.zeros((len(diffusionTimes), nRuns, 2, len(qPoints)))

        for t in range(len(diffusionTimes)):
            diffusionTime = diffusionTimes[t]
            #truePDF = Util.getPDF_sphere(D, diffusionTime, radius, 500, 5)
            #trueCDF = Util.getCDF_sphere(D, diffusionTime, radius, 1000, 5)
            truePDF = Util.getPDF_sphere(D, np.inf, radius, 500, 5)
            trueCDF = Util.getCDF_sphere(D, np.inf, radius, 1000, 5)

            pdfPointsX = np.linspace(0, 1.1*radius, 500)
            pdfPointsY = [truePDF(p) for p in pdfPointsX]

            for r in range(nRuns):
                timeStep = Util.getTimeStep(diffusionTime, nStep)
                l = (6*D*timeStep)**0.5
                envSize = 5*radius
                env = Environment(T2, D, envSize, envSize, envSize)
                part = [Particle3D(*Util.getRandomDirection()*Util.getRandomQuadratic(radius)) for i in range(nPart)]
                #part = [Particle3D(0, 0, 0) for i in range(nPart)]
                sim = Simulation(nStep, timeStep, part, env, [Sphere(0, 0, 0, T2, D, 0, radius)])
                printMessage = "{diffusionTime}ms diffusion time ({t}/{totDt}), run {r}/{nRuns}:"\
                    .format(diffusionTime=diffusionTime, t=t+1, totDt=len(diffusionTimes), r=r+1, nRuns=nRuns)
                histogramTitle = "Probability distribution of the distance travelled by a particle (sphere radius = {radius}um)\n"\
                                 "Diffusion time = {diffusionTime}ms, Number of particles = {nPart}, Number of steps = {n}, run {r}/{nRuns}"\
                    .format(diffusionTime=diffusionTime, nPart=nPart, n=nStep, radius=radius, r=r+1, nRuns=nRuns)
                signals[t, r, :, :] = ValidationTemp.run1SimSignal(sim, printMessage, histogramTitle, trueCDF, pdfPointsX, pdfPointsY, plotHistType, qPoints)
        return signals

    def runPlanesSignalConv(nRuns, plotHistType, diffusionTimes, nStep, nPart, spacing, D, T2, qPoints):
        signals = np.zeros((len(diffusionTimes), nRuns, 2, len(qPoints)))

        for t in range(len(diffusionTimes)):
            diffusionTime = diffusionTimes[t]
            truePDF = Util.getPDF_plane(spacing)
            trueCDF = Util.getCDF_sphere(D, diffusionTime, spacing, 1000, 5)

            pdfPointsX = np.linspace(-0.6*spacing, 0.6*spacing, 500)
            pdfPointsY = [truePDF(p) for p in pdfPointsX]

            for r in range(nRuns):
                timeStep = Util.getTimeStep(diffusionTime, nStep)
                l = (6*D*timeStep)**0.5
                envSize = 4*spacing
                env = Environment(T2, D, envSize, envSize, envSize)
                part = [Particle3D(Util.getRandomU(spacing), Util.getRandomU(envSize), Util.getRandomU(envSize)) for i in range(nPart)]
                #part = [Particle3D(*Util.getRandomDirection()*Util.getRandomQuadratic(radius)) for i in range(nPart)]
                #part = [Particle3D(0, 0, 0) for i in range(nPart)]
                sim = Simulation(nStep, timeStep, part, env, [Planes(T2, D, spacing)])
                printMessage = "{diffusionTime}ms diffusion time ({t}/{totDt}), run {r}/{nRuns}:"\
                    .format(diffusionTime=diffusionTime, t=t+1, totDt=len(diffusionTimes), r=r+1, nRuns=nRuns)
                histogramTitle = "Probability distribution of the distance travelled by a particle (plane spacing = {spacing}um)\n"\
                                 "Diffusion time = {diffusionTime}ms, Number of particles = {nPart}, Number of steps = {n}, run {r}/{nRuns}"\
                    .format(diffusionTime=diffusionTime, nPart=nPart, n=nStep, spacing=spacing, r=r+1, nRuns=nRuns)
                signals[t, r, :, :] = ValidationTemp.run1SimSignal(sim, printMessage, histogramTitle, trueCDF, pdfPointsX, pdfPointsY, plotHistType, qPoints)
        return signals
