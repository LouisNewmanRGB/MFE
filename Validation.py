import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import time
import multiprocessing

from Util import Util
from SimulationResults import SimulationResults
from ValidationCore import ValidationCore

class Validation():

    def runParallel(nRuns, serialFunction, serialArgs, nProcess=8, saveFileName=None):
        rpp = nRuns // nProcess #runs per process
        print("starting {n} runs".format(n=rpp*nProcess))

        manager = multiprocessing.Manager()
        result_dict = manager.dict()
        processes = []
        for number in range(nProcess):
            proc = multiprocessing.Process(target=ValidationCore.runProcess, args=(result_dict, number, rpp, serialFunction, serialArgs,))
            proc.start()
            processes.append(proc)

        for proc in processes:
            proc.join()

        results = np.concatenate(tuple([result_dict[n] for n in range(nProcess)]), axis=-1)

        if saveFileName == None:
            np.save(Util.getFilePath(serialFunction), results)
        else:
            np.save(Util.getFilePath(saveFileName), results)

        return results

    def runValidation(nRuns, diffusionTimes, xDataTuple, xDataType, simulationType, D, T2, parameter=None):
        xData = xDataTuple[0]
        results = np.empty((len(diffusionTimes), len(xData), nRuns), dtype=SimulationResults)
        for t in range(len(diffusionTimes)):
            diffusionTime = diffusionTimes[t]
            for i in range(len(xData)):
                nStep, nPart = ValidationCore.getNStepNPart(xDataTuple, xDataType, i)
                for r in range(nRuns):
                    if xDataType == "comparison":
                        printMessage = "{diffusionTime}ms diffusion time ({t}/{totDt}), {nPart} particles and {n} steps ({i}/{totStep}), run {r}/{nRuns}:"\
                        .format(diffusionTime=diffusionTime, t=t+1, totDt=len(diffusionTimes), nPart=nPart, n=nStep, i=i+1, totStep=len(xData), r=r+1, nRuns=nRuns)
                    elif xDataType == "convergence":
                        printMessage = "{diffusionTime}ms diffusion time ({t}/{totDt}) and {nPart} particles ({i}/{totStep}), run {r}/{nRuns}:"\
                        .format(diffusionTime=diffusionTime, t=t+1, totDt=len(diffusionTimes), nPart=nPart, i=i+1, totStep=len(xData), r=r+1, nRuns=nRuns)
                    else:
                        return print("ERROR: invalid xDataType")

                    timeStep = Util.getTimeStep(diffusionTime, nStep)
                    sim = ValidationCore.getSim(nStep, nPart, timeStep, simulationType, D, T2, parameter)

                    startTime = time.time()
                    sim.run(seed=None, calcData=False)
                    results[t, i, r] = sim.getResults()
                    print(printMessage)
                    print("Computation time: {compTime}s".format(compTime = time.time() - startTime), "\n")
        return results

    def distributionPlot(simResultArray, plotTitle, plotIntermediary, nRuns, diffusionTimes, xDataTuple, xDataType, simulationType, D, parameter=None):
        xData = xDataTuple[0]
        errors = np.zeros((len(diffusionTimes), len(xData), nRuns))
        for t in range(len(diffusionTimes)):
            diffusionTime = diffusionTimes[t]

            if simulationType == "free":
                truePDF = Util.getPDF(D, diffusionTime)
                trueCDF = Util.getCDF(D, diffusionTime)
                pdfPointsX = np.linspace(0, 4 * (6 * D * diffusionTime)**0.5, 500)
                pdfPointsY = truePDF(pdfPointsX)
                histogramTitle = "TODO"
            elif simulationType == "sphere_center":
                radius = parameter
                truePDF = Util.getPDF_sphere(D, diffusionTime, radius, 500, 5)
                trueCDF = Util.getCDF_sphere(D, diffusionTime, radius, 1000, 5)
                pdfPointsX = np.linspace(0, 1.1*radius, 500)
                pdfPointsY = [truePDF(p) for p in pdfPointsX]
                histogramTitle = "TODO"
            else:
                print("ERROR: invalid simulationType")
                return None

            for i in range(len(xData)):
                for r in range(nRuns):
                    simResults = simResultArray[t, i, r]
                    #kolmogorov-smirnov
                    #print(printMessage)
                    test = scipy.stats.kstest(simResults.getDistances(), trueCDF)
                    errors[t, i, r] = test.statistic
                    #print("Kolmogorov-Smirnov test Statistic:", test.statistic)
                    #print("Kolmogorov-Smirnov test pvalue:", test.pvalue, "\n")

                    #histograms
                    if plotIntermediary:
                        bw = 2*scipy.stats.iqr(simResults.getDistances(), rng=(25, 75))/(len(simResults.getDistances()))**(1/3)
                        nBins = int((np.max(simResults.getDistances()) - np.min(simResults.getDistances()))/bw)
                        plt.hist(simResults.getDistances(), bins=nBins, density = True, stacked=True)
                        plt.plot(pdfPointsX, pdfPointsY, color = 'red')
                        plt.legend(["Expected probability density function", "Random walk results histogram"])
                        plt.xlabel("Distance travelled [um]")
                        plt.ylabel("[um-1]")
                        plt.title(histogramTitle)
                        plt.show()

        yLabel = "Supremum distance between exact and empirical CDFs"
        ValidationCore.plotDelegator(errors, diffusionTimes, xData, plotTitle, xDataType, yLabel, pValuePlot=False)

    def meanDisplacementPlot(simResultArray, plotTitle, nRuns, diffusionTimes, xDataTuple, xDataType):
        xData = xDataTuple[0]
        meanDisp = np.zeros((len(diffusionTimes), len(xData), nRuns))
        for t in range(len(diffusionTimes)):
            for i in range(len(xData)):
                for r in range(nRuns):
                    simResults = simResultArray[t, i, r]
                    meanDisp[t, i, r] = np.linalg.norm(simResults.getAvgDisplacements())
        yLabel = "Norm of Mean Displacement"
        ValidationCore.plotDelegator(meanDisp, diffusionTimes, xData, plotTitle, xDataType, yLabel, pValuePlot=True)

    def RMSDplot(simResultArray, plotTitle, nRuns, diffusionTimes, xDataTuple, xDataType, simulationType, D, parameter=None):
        xData = xDataTuple[0]
        errors = np.zeros((len(diffusionTimes), len(xData), nRuns))
        for t in range(len(diffusionTimes)):
            diffusionTime = diffusionTimes[t]

            if simulationType == "sphere_uniform":
                radius = parameter
                expectedMean = Util.RMSD_sphere(radius, D, diffusionTime, 20)
            elif simulationType == "planes":
                spacing = parameter
                expectedMean = Util.RMSD_planes(spacing, D, diffusionTime, 20)
            else:
                return print("ERROR: invalid simulationType")

            for i in range(len(xData)):
                for r in range(nRuns):
                    simResults = simResultArray[t, i, r]
                    if simulationType == "sphere_uniform":
                        errors[t, i, r] = (np.average( simResults.getDistances()**2))**0.5 - expectedMean
                    else:
                        errors[t, i, r] = (np.average( simResults.getDisplacements()[:, 0]**2))**0.5 - expectedMean

        yLabel = "Difference between theoretical and simulated RMSD"
        ValidationCore.plotDelegator(errors, diffusionTimes, xData, plotTitle, xDataType, yLabel, pValuePlot=True)

    def signalPlot(simResultArray, plotIntermediary, plotTitle, qPoints, nRuns, diffusionTimes, xDataTuple, xDataType, simulationType, D, parameter=None, maxStd=None):
        xData = xDataTuple[0]
        RMSEs = np.empty((len(diffusionTimes), len(xData), nRuns))
        qVectors = np.array([[q, 0, 0] for q in qPoints])
        graphTitle = "RWS of diffusion in an impermeable sphere of radius 8um\n" \
                     "diffusion time = 10ms, 38461 particles and 26 time steps (1M total steps)"

        for t in range(len(diffusionTimes)):
            diffusionTime = diffusionTimes[t]
            trueSignalPoints = ValidationCore.getTrueSignalPoints(simulationType, qPoints, D, diffusionTime, parameter)
            for i in range(len(xData)):
                for r in range(nRuns):
                    simResults = simResultArray[t, i, r]

                    simulatedSignal, stdsSample, stdsPopulation = simResults.getSGPSignal(qVectors, includeStd=True)
                    RMSEs[t, i, r] = ( np.average((trueSignalPoints - simulatedSignal)**2) )**0.5

                    if maxStd != None:
                        maxPopStd = np.max(stdsPopulation)
                        nMin = int((maxPopStd/maxStd)**2)
                        print("Number of particles required to reach a sample standard deviation of {maxStd}: {nMin}".format(maxStd=maxStd, nMin=nMin))

                    if plotIntermediary:
                        plt.plot(qPoints, trueSignalPoints, color="r")#colors[t], marker=".")
                        plt.errorbar(qPoints, simulatedSignal, stdsSample, fmt=".", color="g")
                        plt.legend(["Theoretical signal", "Simulated signal", "Simulated signal (real part)"])
                        plt.xlabel("q=gamma*G*delta [um-1]")
                        plt.ylabel("SGP Signal attenuation")
                        plt.title(graphTitle)
                        plt.yscale("log")
                        plt.grid()
                        plt.show()
        yLabel = "RMSE between exact and empirical SGP signals"
        ValidationCore.plotDelegator(RMSEs, diffusionTimes, xData, plotTitle, xDataType, yLabel, pValuePlot=False)

    def averageSignalPlot(simResultArray, plotIntermediary, plotTitle, qPoints, nRuns, diffusionTimes, xDataTuple, xDataType, simulationType, D, parameter=None):
        xData = xDataTuple[0]
        RMSEs = np.empty((len(diffusionTimes), len(xData), 1))
        qVectors = np.array([[q, 0, 0] for q in qPoints])

        for t in range(len(diffusionTimes)):
            diffusionTime = diffusionTimes[t]
            trueSignalPoints = ValidationCore.getTrueSignalPoints(simulationType, qPoints, D, diffusionTime, parameter)

            for i in range(len(xData)):
                signals = np.empty((nRuns, len(qPoints)))
                nStep, nPart = ValidationCore.getNStepNPart(xDataTuple, xDataType, i)
                if simulationType == "sphere_uniform":
                    graphTitle = "{nRuns} run average of SGP signal attenuation in an impermeable sphere (sphere radius = {radius}um)\n"\
                                 "Diffusion time = {diffusionTime}ms, Number of particles = {nPart}, Number of steps = {n}"\
                        .format(radius=parameter, diffusionTime=diffusionTime, nPart=nPart, n=nStep, nRuns=nRuns)
                elif simulationType == "planes":
                    graphTitle = "{nRuns} run average of SGP signal attenuation between two impermeable planes (plane spacing = {spacing}um)\n"\
                                 "Diffusion time = {diffusionTime}ms, Number of particles = {nPart}, Number of steps = {n}"\
                        .format(spacing=parameter, diffusionTime=diffusionTime, nPart=nPart, n=nStep, nRuns=nRuns)
                else:
                    return print("ERROR: invalid simulationType")

                for r in range(nRuns):
                    simResults = simResultArray[t, i, r]
                    signals[r, :] = simResults.getSGPSignal(qVectors)

                signalAverage = np.average(signals, axis=0)
                signalStd = np.std(signals, axis=0)
                RMSEs[t, i, 0] = ( np.average((trueSignalPoints - signalAverage)**2) )**0.5
                #print("Signal root mean square error:", rmse, "\n")

                if plotIntermediary:
                    plt.plot(qPoints, trueSignalPoints, color="r")
                    plt.errorbar(qPoints, signalAverage, signalStd, fmt="o")
                    plt.legend(["Theoretical signal", "Average of simulated signal ({nRuns} runs)".format(nRuns=nRuns)])
                    plt.xlabel("q=gamma*G*delta [um-1]")
                    plt.ylabel("Signal attenuation")
                    plt.title(graphTitle)
                    plt.yscale("log")
                    plt.grid()
                    plt.show()
        yLabel = "RMSE between exact and average empirical SGP signals"
        ValidationCore.plotDelegator(RMSEs, diffusionTimes, xData, plotTitle, xDataType, yLabel, pValuePlot=False)
