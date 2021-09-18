import numpy as np
import random
import math
import cmath
import matplotlib.pyplot as plt
import matplotlib.cm as cm #color map

from Util import Util

class Simulation():
    TIME_TOL = 1e-10 #TODO: improve
    REL_SPACE_TOL = 1e-10

    def __init__(self, nStep, timeStep, particles, environment, compartments=[]):
        self.particles = particles
        self.environment = environment
        self.compartments = compartments + [environment]
        self.nStep = nStep
        self.timeStep = timeStep
        self.startingPositions = np.array([particle.getTruePos() for particle in particles])
        
        #for results
        self.displacements = None
        self.distances = None
        self.signal = None
        
        #for trajectory plotting
        self.positions = None
        self.truePositions = None
        self.posStepIndices = None
        self.truePosStepIndices = None
        
    def getTimeStep(self):
        return self.timeStep
        
    def getPositions(self):
        return self.positions
    
    def getTruePositions(self):
        return self.truePositions
        
    def findCompartment(self, particle):
        for compartment in self.compartments:
            if compartment.contains(particle):
                return compartment
    
    def run(self, seed=None, calcData = False):
        #calcData determines whether intermediate positions will be stored in memory
        
        #reset
        self.diplacements = None
        self.ditances = None
        
        nPart = len(self.particles)
        if seed !=None:
            random.seed(seed)
        if calcData:
            self.positions = [[] for i in range(nPart)]#np.zeros((nPart,self.nStep,3))
            self.truePositions = [[] for i in range(nPart)]#np.zeros((nPart,self.nStep,3))
            self.posStepIndices = [[] for i in range(nPart)] #list of indices of "steps" for every particle
            self.truePosStepIndices = [[] for i in range(nPart)]
        for p in range(nPart):
            particle = self.particles[p]
            particle.setVelocity(Util.getRandomDirection())
            newComp = self.findCompartment(particle)
            particle.changeCompartment(newComp, self.timeStep)
            if calcData:
                self.positions[p].append(particle.getPos().copy())
                self.truePositions[p].append(particle.getTruePos().copy())
                self.posStepIndices[p].append(0)
                self.truePosStepIndices[p].append(0)
            for n in range(self.nStep):
                self.nextStep(p, calcData)
    
    def nextStep(self, particleIndex, calcData):
        particle = self.particles[particleIndex]
        t = 0 #elapsed time during the step
        particle.setVelocity(particle.getSpeed()*Util.getRandomDirection()) #random direction at each step
        
        while t < self.timeStep:
            nComp = len(self.compartments)
            intersections = np.array([None]*nComp)
            reachTimes = np.array([None]*nComp)
            for c in range(nComp):
                intersections[c] = self.compartments[c].findIntersection(particle)
                if type(intersections[c]) == np.ndarray:
                    reachTimes[c] = np.linalg.norm(intersections[c] - particle.getPos())/particle.getSpeed()
                else:
                    reachTimes[c] = math.inf 
            firstIndex = np.argmin(reachTimes)
            reachTime, compartment = reachTimes[firstIndex], self.compartments[firstIndex]

            if t + reachTime < self.timeStep:
                # if there is an intersection
                particle.move(reachTime)
                preCollidePos = particle.getPos().copy()
                compartment.collide(particle, intersections[firstIndex], reachTime, self)
                t += reachTime
                if calcData and (preCollidePos != particle.getPos()).any():
                    self.positions[particleIndex].append(preCollidePos)
            else:
                particle.move(self.timeStep - t)
                t = self.timeStep
            if calcData:
                self.positions[particleIndex].append(particle.getPos().copy())
                self.truePositions[particleIndex].append(particle.getTruePos().copy())
        if calcData:
            self.posStepIndices[particleIndex].append(len(self.positions[particleIndex]) - 1)
            self.truePosStepIndices[particleIndex].append(len(self.truePositions[particleIndex]) - 1)
        
    def plot(self, positionType = True):
        #3D plot of trajectories
        #positionType represents whether to plot "true" or "in cell" positions
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        nPart = len(self.particles)
        colors = cm.rainbow(np.linspace(0, 1, nPart))
        
        if positionType:
            data = self.truePositions
            indices = self.truePosStepIndices
            superMax = Util.recursiveMax(data)
            ax.set_xlim3d([-superMax, superMax])
            ax.set_ylim3d([-superMax, superMax])
            ax.set_zlim3d([-superMax, superMax])
        else:
            ax.set_xlim3d([self.environment.getPos()[0] - self.environment.getSize()[0]/2,
                            self.environment.getPos()[0] + self.environment.getSize()[0]/2])
            ax.set_ylim3d([self.environment.getPos()[1] - self.environment.getSize()[1]/2,
                            self.environment.getPos()[1] + self.environment.getSize()[1]/2])
            ax.set_zlim3d([self.environment.getPos()[2] - self.environment.getSize()[2]/2,
                            self.environment.getPos()[2] + self.environment.getSize()[2]/2])
            data = self.positions
            indices = self.posStepIndices

        for p in range(nPart):
            Util.plotPoints(data[p], ax.plot, colors[p])
            Util.plotPoints([data[p][i] for i in range(len(data[p])) if i in indices[p]], ax.scatter, colors[p])

        for compartment in self.compartments:
            compartment.plot(ax)

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        plt.show()
        
    def getStepLengths(self):
        nPart = len(self.particles)
        stepLengths = [[] for i in range(nPart)]
        for p in range(nPart):
            posList = self.truePositions[p]
            d = 0
            for s in range(1,len(posList)):
                d += np.linalg.norm(posList[s] - posList[s-1])
                if s in self.truePosStepIndices[p]:
                    stepLengths[p].append(d)
                    d = 0
        return stepLengths
        
    def getDistances(self):
        if type(self.distances) != np.array:
            self.distances = np.array([np.linalg.norm(disp) for disp in self.getDisplacements()])
        return self.distances
        
    def getDisplacements(self):
        if type(self.displacements) != np.array:
            endPos = np.array([particle.getTruePos() for particle in self.particles])
            self.displacements = endPos - self.startingPositions
        return self.displacements

    def getSGPSignal(self, gamma, G, delta):
        if self.signal == None:
            nPart = len(self.particles)
            res = 0
            for p in range(nPart):
                res += cmath.exp(1j*gamma*delta*np.dot(self.getDisplacements()[p], G))*self.particles[p].getSignal()
            self.signal = res/nPart
        return self.signal