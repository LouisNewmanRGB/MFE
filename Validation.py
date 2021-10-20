import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy.stats
import time
import multiprocessing

from Particle3D import Particle3D
from Environment import Environment
from Simulation import Simulation
from Sphere import Sphere
from Planes import Planes
from Util import Util


class Validation():

    def runProcess(result_dict, number, rpp, plotHist, serialFunction, serialArgs):
        print("Starting process", number)
        errors = serialFunction(rpp, plotHist, *serialArgs)
        print("Finished process", number)
        result_dict[number] = errors

    def runParallel(nRuns, plotHist, serialFunction, serialArgs, nProcess=8):
        rpp = nRuns // nProcess #runs per process
        print("starting {n} runs".format(n=rpp*nProcess))

        manager = multiprocessing.Manager()
        result_dict = manager.dict()
        processes = []
        for number in range(nProcess):
            proc = multiprocessing.Process(target=Validation.runProcess, args=(result_dict, number, rpp, plotHist, serialFunction, serialArgs,))
            proc.start()
            processes.append(proc)

        for proc in processes:
            proc.join()

        errors = []
        for number in range(nProcess):
            errors.extend(result_dict[number])
        print("errorLength", len(errors))
        return errors
        #np.array(result_dict.values())


    def run1Sim(sim, printMessage, histogramTitle, trueCDF, pdfPointsX, pdfPointsY, plotHist):
        startTime = time.time()
        sim.run(seed=None, calcData=False)
        sim = sim.getResults()
        print("number of particles", len(sim.getDistances()))
        #print("distances", sim.getDistances())

        #kolmogorov-smirnov
        print(printMessage)
        print("Computation time: {compTime}s".format(compTime = time.time() - startTime))
        test = scipy.stats.kstest(sim.getDistances(), trueCDF)
        print("Kolmogorov-Smirnov test Statistic:", test.statistic)
        print("Kolmogorov-Smirnov test pvalue:", test.pvalue, "\n")

        #histograms
        if plotHist:
            bw = 2*scipy.stats.iqr(sim.getDistances(), rng=(25, 75))/(len(sim.getDistances()))**(1/3)
            nBins = int((np.max(sim.getDistances()) - np.min(sim.getDistances()))/bw)
            plt.hist(sim.getDistances(), bins=nBins, density = True, stacked=True)
            plt.plot(pdfPointsX, pdfPointsY, color = 'red')
            plt.legend(["Expected probability density function", "Random walk results histogram"])
            plt.xlabel("Distance travelled [um]")
            plt.ylabel("[um-1]")
            plt.title(histogramTitle)
            plt.show()

        #numbers
        #print("numer:", np.average(np.power(sim.getDistances(), 2)))
        #print("theor:", 6*D*diffusionTime)

        return test.statistic

    def comparisonFinalPlot(errors, diffusionTimes, nStep, plotTitle):
        print("ERRORS:", errors)
        #print("PVALUES:", pvalues)
        averageErrors = np.average(errors, axis=2)
        stdDevs = np.std(errors, axis = 2)
        print("AVERAGE ERRORS:", averageErrors)
        colors = cm.rainbow(np.linspace(0, 1, len(diffusionTimes)))
        for t in range(len(diffusionTimes)):
            #plt.scatter(nStep, averageErrors[t,:], color = colors[t])
            plt.errorbar(nStep, averageErrors[t,:], stdDevs[t,:], fmt="o", color = colors[t], ecolor = colors[t])
        plt.xscale("log")
        plt.xlabel("Number of steps")
        plt.ylabel("Supremum distance between exact and empirical CDFs")
        plt.legend(["Diffusion time = {diffusionTime}ms".format(diffusionTime=dt) for dt in diffusionTimes])
        plt.title(plotTitle)
        plt.show()

    def convergenceFinalPlot(errors, diffusionTimes, nParts, plotTitle):
        print("ERRORS:", errors)
        averageErrors = np.average(errors, axis=2)
        print("AVERAGE ERRORS:", averageErrors)
        stdDevs = np.std(errors, axis = 2)
        for logScale in [False, True]:
            if logScale:
                plt.yscale("log")
                finalPoints = np.array([nParts[0]/10, nParts[-1]*10])
            else:
                finalPoints = 10**np.linspace(np.log10(nParts[0]/1.5), np.log10(nParts[-1]*5), 500)

            colors = cm.rainbow(np.linspace(0, 1, len(diffusionTimes)))
            for t in range(len(diffusionTimes)):
                #plt.scatter(nParts, averageErrors[t,:], color = colors[t])
                plt.errorbar(nParts, averageErrors[t,:], stdDevs[t,:], fmt="o", color = colors[t], ecolor = colors[t])
            plt.plot(finalPoints, finalPoints**(-0.5), color = "black")

            plt.xscale("log")
            plt.xlabel("Number of particles")
            plt.ylabel("Supremum distance between exact and empirical CDFs")
            plt.legend(["Theoretical distances\n(n^-1/2 convergence rate)"] +
                       ["Diffusion time = {diffusionTime}ms".format(diffusionTime=dt) for dt in diffusionTimes])
            plt.title(plotTitle)
            plt.show()

    def runSphereComp(nRuns, plotHist, diffusionTimes, nStep, totalSteps, radius, D, T2):
        errors = np.zeros((len(diffusionTimes), len(nStep), nRuns))
        #pvalues = np.zeros((len(diffusionTimes), len(nStep), nRuns))

        for t in range(len(diffusionTimes)):
            diffusionTime = diffusionTimes[t]
            truePDF = Util.getPDF_sphere(D, diffusionTime, radius, 500, 5)
            trueCDF = Util.getCDF_sphere(D, diffusionTime, radius, 500, 5)
            pdfPointsX = np.linspace(0, 1.1*radius, 500)
            pdfPointsY = [truePDF(p) for p in pdfPointsX] #truePDF(pdfPointsX)
            for i in range(len(nStep)):
                for r in range(nRuns):
                    n = nStep[i]
                    diffusionTime = diffusionTimes[t]
                    timeStep = Util.getTimeStep(diffusionTime, nStep)
                    nPart = int(totalSteps / n)
                    l = (6*D*timeStep)**0.5
                    envSize = 5*radius
                    env = Environment(T2, D, envSize, envSize, envSize)
                    part = [Particle3D(0, 0, 0) for i in range(nPart)]
                    sim = Simulation(n, timeStep, part, env, [Sphere(0, 0, 0, T2, D, 0, radius)])
                    printMessage = "{diffusionTime}ms diffusion time ({t}/{totDt}), {nPart} particles and {n} steps ({i}/{totStep}), run {r}/{nRuns}:"\
                        .format(diffusionTime=diffusionTime, t=t+1, totDt=len(diffusionTimes), nPart=nPart, n=n, i=i+1, totStep=len(nStep), r=r+1, nRuns=nRuns)
                    histogramTitle = "Probability distribution of the distance travelled by a particle (sphere radius = {radius}um)\n"\
                                     "Diffusion time = {diffusionTime}ms, Number of particles = {nPart}, Number of steps = {n} ({totalSteps} particles x steps), run {r}/{nRuns}"\
                        .format(radius=radius, diffusionTime=diffusionTime, nPart=nPart, n=n, totalSteps=totalSteps, r=r+1, nRuns=nRuns)

                    errors[t, i, r] = Validation.run1Sim(sim, printMessage, histogramTitle, trueCDF, pdfPointsX, pdfPointsY, plotHist)
                    #runSim(i, t, r, True)
        return errors

    def runFreeComp(nRuns, plotHist, diffusionTimes, nStep, totalSteps, D, T2):
        errors = np.zeros((len(diffusionTimes), len(nStep), nRuns))

        for t in range(len(diffusionTimes)):
            diffusionTime = diffusionTimes[t]
            truePDF = Util.getPDF(D, diffusionTime)
            trueCDF = Util.getCDF(D, diffusionTime)

            pdfPointsX = np.linspace(0, 4 * (6 * D * diffusionTime)**0.5, 500)
            pdfPointsY = truePDF(pdfPointsX)
            for i in range(len(nStep)):
                for r in range(nRuns):
                    n = nStep[i]
                    diffusionTime = diffusionTimes[t]
                    timeStep = Util.getTimeStep(diffusionTime, nStep)
                    nPart = int(totalSteps / n)
                    l = (6*D*timeStep)**0.5
                    envSize = 10*l
                    env = Environment(T2, D, envSize, envSize, envSize)
                    part = [Particle3D(Util.getRandomU(envSize),Util.getRandomU(envSize),Util.getRandomU(envSize)) for i in range(nPart)]
                    sim = Simulation(n, timeStep, part, env)
                    printMessage = "{diffusionTime}ms diffusion time ({t}/{totDt}), {nPart} particles and {n} steps ({i}/{totStep}), run {r}/{nRuns}:"\
                        .format(diffusionTime=diffusionTime, t=t+1, totDt=len(diffusionTimes), nPart=nPart, n=n, i=i+1, totStep=len(nStep), r=r+1, nRuns=nRuns)
                    histogramTitle = "Probability distribution of the distance travelled by a particle\n"\
                                     "Diffusion time = {diffusionTime}ms, Number of particles = {nPart}, Number of steps = {n}, ({totalSteps} particles x steps), run {r}/{nRuns}"\
                        .format(diffusionTime=diffusionTime, nPart=nPart, n=n, totalSteps=totalSteps, r=r+1, nRuns=nRuns)
                    errors[t, i, r] = Validation.run1Sim(sim, printMessage, histogramTitle, trueCDF, pdfPointsX, pdfPointsY, plotHist)
        return errors


    def runFreeConv(nRuns, plotHist, diffusionTimes, nStep, nParts, D, T2):
        errors = np.zeros((len(diffusionTimes), len(nParts), nRuns))
        pvalues = np.zeros((len(diffusionTimes), len(nParts), nRuns))

        for t in range(len(diffusionTimes)):
            diffusionTime = diffusionTimes[t]
            truePDF = Util.getPDF(D, diffusionTime)
            trueCDF = Util.getCDF(D, diffusionTime)

            pdfPointsX = np.linspace(0, 4 * (6 * D * diffusionTime)**0.5, 500)
            pdfPointsY = truePDF(pdfPointsX)
            for i in range(len(nParts)):
                for r in range(nRuns):
                    nPart = nParts[i]
                    diffusionTime = diffusionTimes[t]
                    timeStep = Util.getTimeStep(diffusionTime, nStep)
                    l = (6*D*timeStep)**0.5
                    envSize = 10*l
                    env = Environment(T2, D, envSize, envSize, envSize)
                    part = [Particle3D(Util.getRandomU(envSize),Util.getRandomU(envSize),Util.getRandomU(envSize)) for i in range(nPart)]
                    sim = Simulation(nStep, timeStep, part, env)
                    printMessage = "{diffusionTime}ms diffusion time ({t}/{totDt}) and {nPart} particles ({i}/{totStep}), run {r}/{nRuns}:"\
                        .format(diffusionTime=diffusionTime, t=t+1, totDt=len(diffusionTimes), nPart=nPart, i=i+1, totStep=len(nParts), r=r+1, nRuns=nRuns)
                    histogramTitle = "Probability distribution of the distance travelled by a particle\n"\
                                     "Diffusion time = {diffusionTime}ms, Number of particles = {nPart}, Number of steps = {n}, run {r}/{nRuns}"\
                        .format(diffusionTime=diffusionTime, nPart=nPart, n=nStep, r=r+1, nRuns=nRuns)
                    errors[t, i, r] = Validation.run1Sim(sim, printMessage, histogramTitle, trueCDF, pdfPointsX, pdfPointsY, plotHist)
        return errors

    def runSphereConv(nRuns, plotHist, diffusionTimes, nStep, nParts, radius, D, T2):
        errors = np.zeros((len(diffusionTimes), len(nParts), nRuns))

        for t in range(len(diffusionTimes)):
            diffusionTime = diffusionTimes[t]
            truePDF = Util.getPDF_sphere(D, diffusionTime, radius, 500, 5)
            trueCDF = Util.getCDF_sphere(D, diffusionTime, radius, 1000, 5)

            pdfPointsX = np.linspace(0, 1.1*radius, 500)
            pdfPointsY = [truePDF(p) for p in pdfPointsX]
            for i in range(len(nParts)):
                for r in range(nRuns):
                    nPart = nParts[i]
                    diffusionTime = diffusionTimes[t]
                    timeStep = Util.getTimeStep(diffusionTime, nStep)
                    l = (6*D*timeStep)**0.5
                    envSize = 5*radius
                    env = Environment(T2, D, envSize, envSize, envSize)
                    #part = [Particle3D(*Util.getRandomDirection()*Util.getRandomQuadratic(radius)) for i in range(nPart)]
                    part = [Particle3D(0, 0, 0) for i in range(nPart)]
                    sim = Simulation(nStep, timeStep, part, env, [Sphere(0, 0, 0, T2, D, 0, radius)])
                    printMessage = "{diffusionTime}ms diffusion time ({t}/{totDt}) and {nPart} particles ({i}/{totStep}), run {r}/{nRuns}:"\
                        .format(diffusionTime=diffusionTime, t=t+1, totDt=len(diffusionTimes), nPart=nPart, i=i+1, totStep=len(nParts), r=r+1, nRuns=nRuns)
                    histogramTitle = "Probability distribution of the distance travelled by a particle (sphere radius = {radius}um)\n"\
                                     "Diffusion time = {diffusionTime}ms, Number of particles = {nPart}, Number of steps = {n}, run {r}/{nRuns}"\
                        .format(diffusionTime=diffusionTime, nPart=nPart, n=nStep, radius=radius, r=r+1, nRuns=nRuns)
                    errors[t, i, r] = Validation.run1Sim(sim, printMessage, histogramTitle, trueCDF, pdfPointsX, pdfPointsY, plotHist)
        return errors

    def run1SimSignal(sim, printMessage, histogramTitle, trueCDF, pdfPointsX, pdfPointsY, plotHistType, qPoints):
        startTime = time.time()
        sim.run(seed=None, calcData=False)
        sim = sim.getResults()
        print("number of particles", len(sim.getDistances()))
        #print("distances", sim.getDistances())

        print(printMessage)
        print("Computation time: {compTime}s".format(compTime = time.time() - startTime))

        if plotHistType != None:
            if plotHistType == "positions_x":
                positions = sim.getEndPositions()[:, 0] #Xs
            elif plotHistType == "positions_norm":
                positions = np.linalg.norm(sim.getEndPositions(), axis=1)
            else:
                print("Error: invalid data requested for histogram")
            bw = 2*scipy.stats.iqr(positions, rng=(25, 75))/(len(positions))**(1/3)
            nBins = int((np.max(positions) - np.min(sim.getDistances()))/bw)
            plt.hist(positions, bins=nBins, density = True, stacked=True)
            plt.plot(pdfPointsX, pdfPointsY, color = 'red')
            plt.legend(["Expected probability density function", "Random walk results histogram"])
            plt.xlabel("Distance travelled [um]")
            plt.ylabel("[um-1]")
            plt.title(histogramTitle)
            plt.show()

        #signal computations:
        simulatedSignal = np.empty(len(qPoints))
        simulatedSignalReal = np.empty(len(qPoints))
        print("Average displacement vector:", sim.getAvgDisplacements())
        for q in range(len(qPoints)):
            qVec = np.array([qPoints[q], 0, 0])
            simulatedSignal[q] = sim.getSGPSignal(qVec, real=False)
            simulatedSignalReal[q] = sim.getSGPSignal(qVec, real=True)
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
                signals[t, r, :, :] = Validation.run1SimSignal(sim, printMessage, histogramTitle, trueCDF, pdfPointsX, pdfPointsY, plotHistType, qPoints)
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
                signals[t, r, :, :] = Validation.run1SimSignal(sim, printMessage, histogramTitle, trueCDF, pdfPointsX, pdfPointsY, plotHistType, qPoints)
        return signals
