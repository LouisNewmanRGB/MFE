import numpy as np
import cmath


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

    def getSGPSignal(self, qVector, real=False):
        #if self.signal == None:
        res = 0
        disp = self.getDisplacements()
        if real:
            for p in range(self.nPart):
                if disp[p][0] > 0:
                    res += 2*np.cos(np.dot(disp[p], qVector))*self.particles[p].getSignal()
        else:
            #for p in range(self.nPart):
            #    res += cmath.exp(1j*np.dot(disp[p], qVector))*self.particles[p].getSignal()
            T2Signal = [p.getSignal() for p in self.particles]
            res = np.dot( np.exp(1j*np.matmul(disp, qVector)), T2Signal)
        return abs(res)/self.nPart
        #return self.signal
