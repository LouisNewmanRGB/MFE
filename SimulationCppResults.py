import numpy as np

class SimulationCppResults():
    #old constructor for non wrapped C++
    """
    def __init__(self, data):
        self.nPart = data.shape[0]
        self.startingPositions = data[:, 0:3]
        self.endPositionsVoxel = data[:, 3:6]
        self.endPositions = data[:, 6:9]
        self.T2Signals = data[:, 9]
        if len(data[0]) == 11:
            self.phases = data[:, 10]
        else:
            self.phases = None

        self.displacements = None
        self.distances = None
    """
    def __init__(self, dataVector, nPart):
        self.nPart = nPart
        numbersPerParticle = int(len(dataVector)/nPart)
        data = np.reshape(dataVector, (nPart, numbersPerParticle))

        self.startingPositions = data[:, 0:3]
        self.endPositionsVoxel = data[:, 3:6]
        self.endPositions = data[:, 6:9]
        self.T2Signals = data[:, 9]
        if numbersPerParticle == 11:
            self.phases = data[:, 10]
        else:
            self.phases = None

        self.displacements = None
        self.distances = None

    def getEndPositions(self):
        return self.endPositions

    def getDisplacements(self):
        if type(self.displacements) != np.array:
            self.displacements = self.endPositions - self.startingPositions
        return self.displacements

    def getDistances(self):
        if type(self.distances) != np.array:
            self.distances = np.array([np.linalg.norm(disp) for disp in self.getDisplacements()])
        return self.distances

    def getAvgDisplacements(self):
        return np.average(self.getDisplacements(), axis=0)

    def getSGPSignal(self, qVectors, includeStd=False):
        disp = self.getDisplacements()
        exponentials = np.exp(1j*np.matmul(qVectors, disp.T)) #replace with cos for real signal
        expTimesT2 = np.multiply([self.T2Signals]*len(qVectors), exponentials)
        signal = np.abs(np.average(expTimesT2, axis=1))
        if includeStd:
            squareMods = np.abs(expTimesT2)**2
            stds = (np.average(squareMods, axis=1) - signal**2)**0.5
            return signal, stds/self.nPart**0.5, stds
        else:
            return signal

    def getFiniteGradientSignalGammaG(self, gammaG, includeStd):
        exponentials = np.multiply(np.exp(1j*gammaG*self.phases), self.T2Signals)
        signal = np.abs(np.average(exponentials, axis=0))
        if includeStd:
            squareMods = np.abs(exponentials)**2
            stds = (np.average(squareMods, axis=0) - signal**2)**0.5
            return signal, stds/self.nPart**0.5, stds
        else:
            return signal

    def getFiniteGradientSignalBValue(self, bValue, Delta, delta, includeStd=False):
        gammaG = (bValue/(Delta - delta/3))**0.5 / delta #only works for Stejskal-Tanner
        #gammaG = self.particles[0].getGradientSequence().getGammaG(bValue)
        return self.getFiniteGradientSignalGammaG(gammaG, includeStd)
