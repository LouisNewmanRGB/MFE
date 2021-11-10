import abc #for abstract classes

class AbstractGradientSequence(metaclass=abc.ABCMeta):
    @abc.abstractmethod
    def __init__(self):
        pass

    def getValue(self, time):
        pass

    def getGammaG(self, bValue):
        pass

    def plot(self):
        pass
