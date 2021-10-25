import numpy as np


class SimulationResults():
    def __init__(self, startingPositions, particles):
        self.startingPositions = startingPositions
        self.particles = particles
        self.nPart = len(self.particles)

        self.displacements = None
        self.endPositions = None
        self.distances = None
        self.signal = None

    def getEndPositions(self):
        if type(self.endPositions) != np.array:
            self.endPositions = np.array([particle.getTruePos() for particle in self.particles])
        return self.endPositions

    def getDistances(self):
        if type(self.distances) != np.array:
            self.distances = np.array([np.linalg.norm(disp) for disp in self.getDisplacements()])
        return self.distances

    def getDisplacements(self):
        if type(self.displacements) != np.array:
            #endPos = np.array([particle.getTruePos() for particle in self.particles])
            self.displacements = self.getEndPositions() - self.startingPositions
        return self.displacements

    def getAvgDisplacements(self):
        return np.average(self.getDisplacements(), axis=0)

    def getSGPSignalSlow(self, qVector, real=False):
        """"slower version"""
        #if self.signal == None:
        disp = self.getDisplacements()
        T2Signal = np.array([p.getSignal() for p in self.particles])
        if real:
            qNorm = np.linalg.norm(qVector)
            dispX = disp[:, 0]
            temp = np.multiply(np.cos(qNorm*dispX), T2Signal) #Xs
            res = 2*np.sum(np.where(dispX > 0, temp, 0))
            #for p in range(self.nPart):
            #    if disp[p][0] > 0:
            #        res += 2*np.cos(np.dot(disp[p], qVector))*self.particles[p].getSignal()
        else:
            #for p in range(self.nPart):
            #    res += cmath.exp(1j*np.dot(disp[p], qVector))*self.particles[p].getSignal()
            res = np.dot( np.exp(1j*np.matmul(disp, qVector)), T2Signal)
            #res = res.imag
        return abs(res)/self.nPart
        #return self.signal

    def getSGPSignal(self, qVectors, real=False):
        #if self.signal == None:
        disp = self.getDisplacements()
        T2Signal = np.array([p.getSignal() for p in self.particles])
        if real:
            exponentials = np.cos(np.matmul(qVectors, disp.T))
        else:
            exponentials = np.exp(1j*np.matmul(qVectors, disp.T))

        expTimesT2 = np.multiply([T2Signal]*len(qVectors), exponentials)
        return np.abs(np.average(expTimesT2, axis=1))
        #return self.signal
