import abc #for abstract classes

class AbstractGradientSequence(metaclass=abc.ABCMeta):
    @abc.abstractmethod
    def __init__(self):
        pass

    def getValue(self, time):
        pass

    def getPhaseStep(self, particle, dt, t0, discreteSequenceIndex):
        pass

    def getGammaG(self, bValue):
        pass
