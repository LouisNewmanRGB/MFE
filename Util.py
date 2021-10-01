import numpy as np
import random
import math
import scipy.stats

class Util():

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

    def genAlphaList(nTerms, nIter, a):
        """positive roots (solving for alpha) of a*alpha = tan(a*alpha)"""
        alphaList = [(2*k + 1)*np.pi/2 for k in range(nTerms)]
        for j in range(nIter):
            for k in range(nTerms):
                alphaList[k] = np.arctan(alphaList[k]) + k*np.pi
        return (np.array(alphaList)/a)[1:]

    def getPDF_sphere(D, diffusionTime, radius, nTerms, nIter):
        alphaList = Util.genAlphaList(nTerms, nIter, radius)

        def f(r):
            if r < radius:
                return 3*r**2/(radius**3) + \
                   2*r*np.sum([np.exp(-D*diffusionTime*alpha**2)*np.sin(alpha*r)*alpha/(np.sin(alpha*radius)**2) for alpha in alphaList])/radius
            return 0
        return f

    def getCDF_sphere(D, diffusionTime, radius, nTerms, nIter):
        alphaList = Util.genAlphaList(nTerms, nIter, radius)

        def F(r):

            return np.heaviside(r - radius, 0) + np.heaviside(radius - r, 1)*( (r/radius)**3 + \
                        2*np.sum([np.exp(-D*diffusionTime*alpha**2)* \
                                  (np.sin(alpha*r) - alpha*r*np.cos(alpha*r))/(alpha*np.sin(alpha*radius)**2) for alpha in alphaList], axis=0)/radius )
        return F
