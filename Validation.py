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
from Util import Util

class Validation():

    def runProcess(result_dict, number, rpp, serialFunction, serialArgs):
        print("Starting process", number)
        results = serialFunction(rpp, *serialArgs)
        print("Finished process", number)
        result_dict[number] = results

    def runParallel(nRuns, serialFunction, serialArgs, nProcess=8):
        rpp = nRuns // nProcess #runs per process
        print("starting {n} runs".format(n=rpp*nProcess))

        manager = multiprocessing.Manager()
        result_dict = manager.dict()
        processes = []
        for number in range(nProcess):
            proc = multiprocessing.Process(target=Validation.runProcess, args=(result_dict, number, rpp, serialFunction, serialArgs,))
            proc.start()
            processes.append(proc)

        for proc in processes:
            proc.join()

        results = np.concatenate(tuple([result_dict[n] for n in range(nProcess)]), axis=-1)
        np.save(Util.getFilePath(serialFunction), results)

        return results


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

        return test.statistic

    def comparisonPlot(data, diffusionTimes, nStep, plotTitle, yLabel):
        averageData = np.average(data, axis=2)
        stdDevs = np.std(data, axis = 2)
        colors = cm.rainbow(np.linspace(0, 1, len(diffusionTimes)))
        fig, ax = plt.subplots()
        for t in range(len(diffusionTimes)):
            ax.errorbar(nStep, averageData[t,:], stdDevs[t,:], fmt="o", color = colors[t], ecolor = colors[t])
        plt.xscale("log")
        if np.all(data > 0):
            ax.set_ylim(ymin=0)
        plt.xlabel("Number of steps")
        plt.ylabel(yLabel)
        plt.legend(["Diffusion time = {diffusionTime}ms".format(diffusionTime=dt) for dt in diffusionTimes])
        plt.title(plotTitle)
        plt.grid()
        plt.show()

    def compPValuePlot(pValues, diffusionTimes, nStep, plotTitle):
        colors = cm.rainbow(np.linspace(0, 1, len(diffusionTimes)))
        fig, ax = plt.subplots()
        for t in range(len(diffusionTimes)):
            ax.scatter(nStep, pValues[t,:], color = colors[t])
        plt.xscale("log")
        ax.set_ylim(ymin=0)
        plt.xlabel("Number of steps")
        plt.ylabel("p-value")
        plt.legend(["Diffusion time = {diffusionTime}ms".format(diffusionTime=dt) for dt in diffusionTimes])
        plt.title(plotTitle)
        plt.grid()
        plt.show()

    def meanTestCompPlot(data, diffusionTimes, nStep, plotTitle, yLabel):
        Validation.comparisonPlot(data, diffusionTimes, nStep, plotTitle, yLabel)
        pValues = np.empty((len(diffusionTimes), len(nStep)))
        for t in range(len(diffusionTimes)):
            for i in range(len(nStep)):
                pValues[t, i] = scipy.stats.ttest_1samp(data[t, i, :], 0).pvalue
        plotTitleExtra = "\n" + yLabel + ": Student Test for Null Mean"
        Validation.compPValuePlot(pValues, diffusionTimes, nStep, plotTitle + plotTitleExtra)

    def meanTestCompMD(meanDisps, diffusionTimes, nStep, plotTitle):
        Validation.meanTestCompPlot(meanDisps, diffusionTimes, nStep, plotTitle, "Norm of Mean Displacement")

    def meanTestCompRMSD(RMSDs, diffusionTimes, plotTitle, nStep, radius, D):
        data = RMSDs.copy()
        for t in range(len(diffusionTimes)):
            expectedMean = Util.RMSD_sphere(radius, D, diffusionTimes[t], 20, 5)
            data[t, :, :] = data[t, :, :] - expectedMean
        Validation.meanTestCompPlot(data, diffusionTimes, nStep, plotTitle, "Difference between theoretical and simulated RMSD")

    def diffusionCompPlot(errors, diffusionTimes, nStep, plotTitle):
        Validation.comparisonPlot(errors, diffusionTimes, nStep, plotTitle, yLabel="Supremum distance between exact and empirical CDFs")

    def signalCompPlot(errors, diffusionTimes, nStep, plotTitle):
        Validation.comparisonPlot(errors, diffusionTimes, nStep, plotTitle, yLabel="RMSE between exact and empirical SGP signals")

    def convergenceFinalPlot(errors, diffusionTimes, nParts, plotTitle):
        averageErrors = np.average(errors, axis=2)
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
                    timeStep = Util.getTimeStep(diffusionTime, n)
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
                    timeStep = Util.getTimeStep(diffusionTime, n)
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

    def run1SimSignal(sim, printMessage, qPoints, trueSignalPoints, graphTitle, plotGraph):
        startTime = time.time()
        sim.run(seed=None, calcData=False)
        sim = sim.getResults()
        print(printMessage)
        print("Computation time: {compTime}s".format(compTime = time.time() - startTime))

        #distance computations:
        md = np.linalg.norm( sim.getAvgDisplacements())
        rmsd = (np.average( sim.getDistances()**2))**0.5
        print("Mean displacement:", md)
        print("Root mean square displacement:", rmsd)

        #signal computations:
        qVectors = np.array([[q, 0, 0] for q in qPoints])
        simulatedSignal = sim.getSGPSignal(qVectors, real=False)
        rmse = ( np.average((trueSignalPoints - simulatedSignal)**2) )**0.5
        print("Signal root mean square error:", rmse, "\n")

        if plotGraph:
            simulatedSignalReal = sim.getSGPSignal(qVectors, real=True)
            plt.plot(qPoints, trueSignalPoints, color="r")#colors[t], marker=".")
            plt.plot(qPoints, simulatedSignal, color="g")#colors[t], marker="o")
            plt.plot(qPoints, simulatedSignalReal, color="b")#colors[t], marker="v")
            plt.legend(["Theoretical signal", "Simulated signal", "Simulated signal (real part)"])
            plt.xlabel("q=gamma*G*delta [um-1]")
            plt.ylabel("Signal attenuation")
            plt.title(graphTitle)
            plt.yscale("log")
            plt.grid()
            plt.show()

        return np.array([rmse, md, rmsd])

    def runSphereSignalComp(nRuns, plotGraph, diffusionTimes, nStep, totalSteps, radius, D, T2, qPoints):
        results = np.zeros((len(diffusionTimes), len(nStep), 3, nRuns))

        for t in range(len(diffusionTimes)):
            diffusionTime = diffusionTimes[t]

            trueSignal = Util.getSignal_sphere_fin(radius, D, diffusionTime, 20, 20)
            trueSignalPoints = trueSignal(qPoints)

            for i in range(len(nStep)):
                for r in range(nRuns):
                    n = nStep[i]
                    nPart = int(totalSteps / n)
                    timeStep = Util.getTimeStep(diffusionTime, n)
                    l = (6*D*timeStep)**0.5
                    envSize = 5*radius
                    env = Environment(T2, D, envSize, envSize, envSize)
                    part = [Particle3D(*Util.getRandomDirection()*Util.getRandomQuadratic(radius)) for i in range(nPart)]
                    sim = Simulation(n, timeStep, part, env, [Sphere(0, 0, 0, T2, D, 0, radius)])
                    printMessage = "{diffusionTime}ms diffusion time ({t}/{totDt}), {nPart} particles and {n} steps ({i}/{totStep}), run {r}/{nRuns}:"\
                        .format(diffusionTime=diffusionTime, t=t+1, totDt=len(diffusionTimes), nPart=nPart, n=n, i=i+1, totStep=len(nStep), r=r+1, nRuns=nRuns)
                    graphTitle = "SGP signal attenuation in an impermeable sphere (sphere radius = {radius}um)\n"\
                                     "Diffusion time = {diffusionTime}ms, Number of particles = {nPart}, Number of steps = {n}, run {r}/{nRuns}"\
                        .format(diffusionTime=diffusionTime, nPart=nPart, n=n+1, radius=radius, r=r+1, nRuns=nRuns)
                    results[t, i, :, r] = Validation.run1SimSignal(sim, printMessage, qPoints, trueSignalPoints, graphTitle, plotGraph)
        return results
