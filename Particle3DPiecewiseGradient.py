import numpy as np
import math

from Particle3D import Particle3D

class Particle3DPiecewiseGradient(Particle3D):

    def __init__(self, x, y, z, gradientSequence):
        super(Particle3DPiecewiseGradient, self).__init__(float(x), float(y), float(z))
        self.gradientSequence = gradientSequence
        self.phase = 0
        self.criticalTimes = gradientSequence.getCriticalTimes().copy()
        self.criticalTime = self.criticalTimes.pop(-1)
        self.gradientValue = self.gradientSequence.getCriticalValue(self.criticalTime)
        self.totalElapsedTime = 0

    def move(self, time):
        t = self.totalElapsedTime
        tFinal = t + time
        while len(self.criticalTimes) > 0 and self.criticalTime < tFinal:
            dt = self.criticalTime - t
            #print(dt, self.gradientValue)
            if np.any(self.gradientValue != 0):
                self.phase += self.integrateGradient(dt)
            t = self.criticalTime
            self.criticalTime = self.criticalTimes.pop(-1)
            self.gradientValue = self.gradientSequence.getCriticalValue(self.criticalTime)
            self.truePos += dt*self.velocity
        dtFinal = tFinal - t
        #print(dtFinal, self.gradientValue)
        if np.any(self.gradientValue != 0):
            self.phase += self.integrateGradient(dtFinal)
        self.pos += time*self.velocity
        self.truePos += dtFinal*self.velocity
        self.T2Signal *= math.exp(-time/self.compartment.getT2()) #T2 signal attenuation
        self.totalElapsedTime += time

    def getSignal(self, gammaG):
        return np.exp(1j*gammaG*self.phase)*self.getT2Signal()

    def getGradientSequence(self):
        return self.gradientSequence

    def integrateGradient(self, time):
        return time*np.dot(self.gradientValue, self.truePos + self.velocity*time/2)
