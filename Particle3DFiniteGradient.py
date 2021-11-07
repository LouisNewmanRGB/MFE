import numpy as np
import math

from Particle3D import Particle3D

class Particle3DFiniteGradient(Particle3D):

    def __init__(self, x, y, z, gradientSequence):
        super(Particle3DFiniteGradient, self).__init__(float(x), float(y), float(z))
        self.gradientSequence = gradientSequence
        self.phase = 0
        self.discreteSequenceIndex = 0 #index of a "critical" point in a piecewize gradientSequence
        self.totalElapsedTime = 0

    def move(self, time):
        self.pos += time*self.velocity
        self.truePos += time*self.velocity
        self.T2Signal *= math.exp(-time/self.compartment.getT2()) #T2 signal attenuation
        self.phase += self.gradientSequence.getPhaseStep(self, time, self.totalElapsedTime, self.discreteSequenceIndex)
        self.totalElapsedTime += time

    def incrementDiscreteSequenceIndex(self):
        self.discreteSequenceIndex += 1

    def getSignal(self, gammaG):
        return np.exp(1j*gammaG*self.phase)*self.getT2Signal()

    def getGradientSequence(self):
        return self.gradientSequence
