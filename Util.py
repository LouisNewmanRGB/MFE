import numpy as np
import random
import math
import scipy.stats

class Util():
    TOL = 1e-10

    def rootsReal(array):
        roots = np.roots(array)
        return roots.real[abs(roots.imag)<Util.TOL]

    def recursiveMax(listOfLists):
        if type(listOfLists[0]) == list:
            return np.max([Util.recursiveMax(e) for e in listOfLists])
        else:
            return(np.max(np.abs(listOfLists)))

    def getRandomDirection():
        theta = random.random() * 2*math.pi
        phi = random.random() * math.pi
        return np.array([math.sin(phi) * math.cos(theta), math.sin(phi) * math.sin(theta), math.cos(phi)])

    def getRandomU(size):
        return (random.random()-0.5)*size

    def plotPoints(data, plotter, color):
        xArray = np.array([pos[0] for pos in data])
        yArray = np.array([pos[1] for pos in data])
        zArray = np.array([pos[2] for pos in data])
        plotter(xArray, yArray, zArray, color = color)

    def getPDF(D, diffusionTime):
        def f(x):
            return scipy.stats.maxwell.pdf(x, loc = 0, scale = (2*D*diffusionTime)**0.5)
        return f

    def getCDF(D, diffusionTime):
        def F(x):
            return scipy.stats.maxwell.cdf(x, loc = 0, scale = (2*D*diffusionTime)**0.5)
        return F
