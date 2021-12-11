import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy.stats

from SimulationCpp import Simulation
from Util import Util

class ValidationCore():

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

    def pValuePlot(pValues, diffusionTimes, nStep, plotTitle, xLabel):
        colors = cm.rainbow(np.linspace(0, 1, len(diffusionTimes)))
        fig, ax = plt.subplots()
        for t in range(len(diffusionTimes)):
            ax.scatter(nStep, pValues[t,:], color = colors[t])
        plt.xscale("log")
        ax.set_ylim(ymin=0, ymax=1)
        plt.xlabel(xLabel)
        plt.ylabel("p-value")
        plt.legend(["Diffusion time = {diffusionTime}ms".format(diffusionTime=dt) for dt in diffusionTimes])
        plt.title(plotTitle)
        plt.grid()
        plt.show()

    def convergencePlot(errors, diffusionTimes, nParts, plotTitle, yLabel):
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
            plt.ylabel(yLabel)
            plt.legend(["Theoretical distances\n(n^-1/2 convergence rate)"] +
                       ["Diffusion time = {diffusionTime}ms".format(diffusionTime=dt) for dt in diffusionTimes])
            plt.title(plotTitle)
            plt.grid()
            plt.show()

    def plotDelegator(data, diffusionTimes, xData, plotTitle, xDataType, yLabel, pValuePlot):
        if xDataType == "comparison":
            xLabel = "Number of Steps"
            ValidationCore.comparisonPlot(data, diffusionTimes, xData, plotTitle, yLabel)

        elif xDataType == "convergence":
            xLabel = "Number of Particles"
            ValidationCore.convergencePlot(data, diffusionTimes, xData, plotTitle, yLabel)

        else:
            return print("ERROR: invalid xDataType")

        if pValuePlot:
            pValues = np.empty((len(diffusionTimes), len(xData)))
            for t in range(len(diffusionTimes)):
                for i in range(len(xData)):
                    pValues[t, i] = scipy.stats.ttest_1samp(data[t, i, :], 0).pvalue
            plotTitleExtra = "\n" + yLabel + ": Student Test for Null Mean"
            ValidationCore.pValuePlot(pValues, diffusionTimes, xData, plotTitle + plotTitleExtra, xLabel)

    def getNStepNPart(xDataTuple, xDataType, i):
        if xDataType == "comparison":
            nSteps, totalSteps = xDataTuple[0], xDataTuple[1]
            nStep = nSteps[i]
            nPart = int(totalSteps / nStep)
        elif xDataType == "convergence":
            nParts, nStep = xDataTuple[0], xDataTuple[1]
            nPart = nParts[i]
        else:
            return print("ERROR: invalid xDataType")
        return nStep, nPart

    def getSim(nStep, nPart, timeStep, simulationType, D, T2, parameter=None):
        if simulationType == "free":
            l = (6*D*timeStep)**0.5
            envSize = 10*l
            sim = Simulation(nStep, timeStep, T2, D, envSize, envSize, envSize)
            sim.createStartingPositions(nPart, False)
            #sim.createStartingPositionsAtPos(nPart, 0, 0, 0)
            """
            part = [Particle3D(np.random.uniform(-envSize/2, envSize/2), np.random.uniform(-envSize/2, envSize/2), \
                               np.random.uniform(-envSize/2, envSize/2)) for i in range(nPart)]
            #part = [Particle3D(Util.getRandomU(envSize),Util.getRandomU(envSize),Util.getRandomU(envSize)) for i in range(nPart)]
            compartments = []
            """

        elif simulationType == "sphere_center":
            radius = parameter
            envSize = 5*radius
            sim = Simulation(nStep, timeStep, T2, D, envSize, envSize, envSize)
            sim.addSphere(0, 0, 0, T2, D, 0, radius)
            sim.createStartingPositionsAtPos(nPart, 0, 0, 0)
            """
            part = [Particle3D(0, 0, 0) for i in range(nPart)]
            compartments = [Sphere(0, 0, 0, T2, D, 0, radius)]
            """

        elif simulationType == "sphere_uniform":
            radius = parameter
            envSize = 5*radius
            sim = Simulation(nStep, timeStep, T2, D, envSize, envSize, envSize)
            sim.addSphere(0, 0, 0, T2, D, 0, radius)
            sim.createStartingPositions(nPart, True)
            """
            part = [Particle3D(*Util.getRandomDirection()*Util.getRandomQuadratic(radius)) for i in range(nPart)]
            compartments = [Sphere(0, 0, 0, T2, D, 0, radius)]
            """

        elif simulationType == "planes":
            spacing = parameter
            envSize = 4*spacing
            sim = Simulation(nStep, timeStep, T2, D, envSize, envSize, envSize)
            sim.addPlanes(T2, D, spacing)
            sim.createStartingPositions(nPart, True)
            """
            part = [Particle3D(np.random.uniform(-spacing/2, spacing/2), np.random.uniform(-envSize/2, envSize/2), \
                               np.random.uniform(-envSize/2, envSize/2)) for i in range(nPart)]
            #part = [Particle3D(Util.getRandomU(spacing), Util.getRandomU(envSize), Util.getRandomU(envSize)) for i in range(nPart)]
            compartments = [Planes(T2, D, spacing)]
            """

        else:
            return print("ERROR: invalid simulationType")

        """
        env = Environment(T2, D, envSize, envSize, envSize)
        return Simulation(nStep, timeStep, part, env, compartments)
        """
        return sim

    def getTrueSignalPoints(simulationType, qPoints, D, diffusionTime, parameter):
        if simulationType == "sphere_uniform":
            radius = parameter
            trueSignal = Util.getSignal_sphere_fin(radius, D, diffusionTime, 20, 20, 10)
        elif simulationType == "planes":
            spacing = parameter
            trueSignal = Util.getSignal_plane_fin(spacing, D, diffusionTime, 20)
        else:
            return print("ERROR: invalid simulationType")
        return trueSignal(qPoints)

    """
    def runProcess(result_dict, number, rpp, serialFunction, serialArgs):
        print("Starting process", number)
        results = serialFunction(rpp, *serialArgs)
        print("Finished process", number)
        result_dict[number] = results
    """
