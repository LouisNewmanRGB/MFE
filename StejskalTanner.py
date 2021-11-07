import numpy as np

from AbstractGradientSequence import AbstractGradientSequence

class StejskalTanner(AbstractGradientSequence):
    def __init__(self, t, TE, delta, Delta):
        self.t = t
        self.TE = TE
        self.delta = delta
        self.Delta = Delta
        assert(TE/2 >= t + delta and TE/2 <= t + Delta) #insures that the 180deg pulse is between the gradients
        self.criticalTimes = np.array(list(set([t, t + delta, t + Delta, t + delta + Delta, TE])))
        self.criticalTimes.sort()
        self.criticalValues = np.array([self.getValue(tim) for tim in self.criticalTimes])
        self.maxIndex = len(self.criticalValues) - 1
        print(self.criticalTimes)

    def getValue(self, time):
        value = np.heaviside(time - self.t, 0) - np.heaviside(time - self.t - self.delta, 0) \
                       - np.heaviside(time - self.t - self.Delta, 0) + np.heaviside(time - self.t - self.Delta -self.delta, 0)
        return np.array([value, 0, 0])

    def getPhaseStep(self, particle, dt, t0, discreteSequenceIndex):
        t = t0
        tFinal = t0 + dt
        phase = 0
        criticalIndex = discreteSequenceIndex
        while criticalIndex < self.maxIndex and self.criticalTimes[criticalIndex] < tFinal:
            phase += self.integrateGradient(particle, self.criticalTimes[criticalIndex] - t, self.criticalValues[criticalIndex])
            t = self.criticalTimes[criticalIndex]
            criticalIndex += 1
            particle.incrementDiscreteSequenceIndex()
        phase += self.integrateGradient(particle, tFinal - t, self.criticalValues[criticalIndex])
        return phase

    def getGammaG(self, bValue):
        return (bValue/(self.Delta - self.delta/3))**0.5 / self.delta

    def integrateGradient(self, particle, time, Gvalue):
        #print(time*Gvalue)
        if np.any(Gvalue != 0):
            return time*np.dot(Gvalue, particle.getPos() + particle.getVelocity()*time/2)
        else:
            return 0
