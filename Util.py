import numpy as np
import random
import math

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

    def plotPoints(data, plotter, color):
        xArray = np.array([pos[0] for pos in data])
        yArray = np.array([pos[1] for pos in data])
        zArray = np.array([pos[2] for pos in data])
        plotter(xArray, yArray, zArray, color = color)
