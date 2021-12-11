import numpy as np
import matplotlib.pyplot as plt

from SimulationPython.Simulation.AbstractGradientSequence import AbstractGradientSequence

class StejskalTanner(AbstractGradientSequence):
    def __init__(self, t, TE, delta, Delta):
        self.t = t
        self.TE = TE
        self.delta = delta
        self.Delta = Delta
        assert(TE/2 >= t + delta and TE/2 <= t + Delta) #insures that the 180deg pulse is between the gradients
        self.criticalTimes = list(set([t, t + delta, t + Delta, t + delta + Delta, TE]))
        self.criticalTimes.sort(reverse=True)
        if t==0:
            self.criticalTimes.pop(-1)
        self.criticalValues = {tim:self.getValue(tim) for tim in self.criticalTimes}
        self.maxIndex = len(self.criticalValues) - 1

    def getValue(self, time):
        value = - np.heaviside(time - self.t, 0) + np.heaviside(time - self.t - self.delta, 0) \
                       + np.heaviside(time - self.t - self.Delta, 0) - np.heaviside(time - self.t - self.Delta -self.delta, 0)
        return np.array([value, 0, 0])

    def getGammaG(self, bValue):
        return (bValue/(self.Delta - self.delta/3))**0.5 / self.delta

    def plot(self):
        xMaxs = np.flip(self.criticalTimes)
        xMins = [0] + list(xMaxs[:-1])
        ys = [self.getValue(tim)[0] for tim in xMaxs]
        plt.hlines(ys, xmin=xMins, xmax=xMaxs)
        plt.scatter(xMaxs, ys)
        plt.xlabel("time [ms]")
        plt.ylabel("x Component of magnetic field gradient [au]")
        plt.title("Stejskal-Tanner sequence of parameters \n"
                  "t={t}ms, TE={TE}ms, Delta={Delta}ms, delta={delta}ms".format(t=self.t, TE=self.TE, Delta=self.Delta, delta=self.delta))
        plt.grid()
        plt.show()

    def getCriticalTimes(self):
        return self.criticalTimes

    def getCriticalValue(self, criticalTime):
        return self.criticalValues[criticalTime]

