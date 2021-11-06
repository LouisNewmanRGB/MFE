import numpy as np
import random
import matplotlib.pyplot as plt
import matplotlib.cm as cm #color map
import time

from Util import Util
from SimulationResults import SimulationResults

class Simulation():
    TOL = 1e-10 #relative tolerance

    def __init__(self, nStep, timeStep, particles, environment, compartments=[]):
        self.particles = particles
        self.nPart = len(self.particles)
        self.environment = environment
        self.compartments = compartments + [environment]
        self.nComp = len(self.compartments)
        self.nStep = nStep
        self.timeStep = timeStep
        self.timeTolerance = self.timeStep*Simulation.TOL #time used for particles to move right after collision
        self.startingPositions = np.array([particle.getTruePos() for particle in particles])
        
        #for results
        self.results = None
        
        #for trajectory plotting
        self.positions = None
        self.truePositions = None
        self.dataCompartments = None
        self.posStepIndices = None
        self.truePosStepIndices = None

    def getParticles(self):
        return self.particles

    def getResults(self):
        return self.results
        
    def getTimeStep(self):
        return self.timeStep
        
    def getPositions(self):
        return self.positions
    
    def getTruePositions(self):
        return self.truePositions
        
    def findCompartment(self, pos, excludedComp = None):
        for compartment in self.compartments:
            if compartment != excludedComp and compartment.contains(pos):
                return compartment
    
    def run(self, seed=None, calcData = False, partPrintNumber = None):
        #calcData determines whether intermediate positions will be stored in memory
        for comp in self.compartments[:-1]:
            parent = self.findCompartment(comp.getPos(), excludedComp = comp)
            comp.setParentComp(parent)

        if calcData:
            self.positions = [[] for i in range(self.nPart)]
            self.truePositions = [[] for i in range(self.nPart)]
            self.dataCompartments = [[] for i in range(self.nPart)]
            self.posStepIndices = [[] for i in range(self.nPart)]
            self.truePosStepIndices = [[] for i in range(self.nPart)]

        if seed !=None:
            random.seed(seed)

        startTime = time.time()
        for p in range(self.nPart):
            particle = self.particles[p]
            particle.setVelocity(Util.getRandomDirection())
            newComp = self.findCompartment(particle.getPos())
            particle.changeCompartment(newComp, self.timeStep)
            if calcData:
                self.positions[p].append(particle.getPos().copy())
                self.truePositions[p].append(particle.getTruePos().copy())
                self.dataCompartments[p].append(particle.getCompartment())
                self.posStepIndices[p].append(0)
                self.truePosStepIndices[p].append(0)
            for n in range(self.nStep):
                self.nextStep(p, calcData)
            if partPrintNumber != None and (p+1)%partPrintNumber == 0:
                print("Particle {p}/{nPart}\n{time}s\n".format(p=p+1, nPart=self.nPart, time=time.time() - startTime))
        self.results = SimulationResults(self.startingPositions, self.particles)
    
    def nextStep(self, particleIndex, calcData):
        particle = self.particles[particleIndex]
        t = 0 #elapsed time during the step
        particle.setVelocity(particle.getSpeed()*Util.getRandomDirection()) #random direction at each step
        while t < self.timeStep:
            ray = np.array([particle.getPos(), particle.getVelocity()/particle.getSpeed()])
            distances = np.empty(self.nComp)
            maxDistance = particle.getSpeed()*(self.timeStep - t)
            for c in range(self.nComp):
                distances[c] = self.compartments[c].findIntersection(ray, maxDistance)
            firstIndex = np.argmin(distances)
            compartment = self.compartments[firstIndex]
            reachTime = distances[firstIndex]/particle.getSpeed()
            
            if t + reachTime + self.timeTolerance < self.timeStep:
                # if there is an intersection
                oldPos = particle.getPos().copy()
                particle.move(reachTime)
                preCollidePos = particle.getPos().copy()
                compartment.collide(particle, oldPos, ray[0] + distances[firstIndex]*ray[1], self)
                if calcData and (preCollidePos != particle.getPos()).any():
                    self.positions[particleIndex].append(preCollidePos)
                particle.move(self.timeTolerance) #moving to avoid ray origin intersecting with compartments on next iteration
                t += reachTime + self.timeTolerance
            else:
                particle.move(self.timeStep - t)
                t = self.timeStep
            if calcData:
                self.positions[particleIndex].append(particle.getPos().copy())
                self.truePositions[particleIndex].append(particle.getTruePos().copy())
        if calcData:
            self.dataCompartments[particleIndex].append(particle.getCompartment())
            self.posStepIndices[particleIndex].append(len(self.positions[particleIndex]) - 1)
            self.truePosStepIndices[particleIndex].append(len(self.truePositions[particleIndex]) - 1)

    def plot(self, positionType = True):
        #3D plot of trajectories
        #positionType represents whether to plot "true" or "in cell" positions
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        colors = cm.rainbow(np.linspace(0, 1, self.nPart))
        compartmentColors = cm.gnuplot(np.linspace(0, 1, self.nComp))
        
        if positionType:
            data = self.truePositions
            indices = self.truePosStepIndices
            superMax = Util.recursiveMax(data)
            ax.set_xlim3d([-superMax, superMax])
            ax.set_ylim3d([-superMax, superMax])
            ax.set_zlim3d([-superMax, superMax])
        else:
            ax.set_xlim3d([self.environment.getSize()[0]/2, -self.environment.getSize()[0]/2])
            ax.set_ylim3d([self.environment.getSize()[1]/2, -self.environment.getSize()[1]/2])
            ax.set_zlim3d([self.environment.getSize()[2]/2, -self.environment.getSize()[2]/2])
            data = self.positions
            indices = self.posStepIndices

        for p in range(self.nPart):
            Util.plotPoints(data[p], ax.plot, colors[p]) #line plot

            #dot plot
            plotData = [data[p][i] for i in range(len(data[p])) if i in indices[p]]
            for c in range(self.nComp):
                Util.plotPoints([plotData[j] for j in range(len(plotData)) \
                                 if self.dataCompartments[p][j] == self.compartments[c]], ax.scatter, compartmentColors[c])

        for compartment in self.compartments:
            compartment.plot(ax)

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        plt.show()
        
    def getStepLengths(self):
        stepLengths = [[] for i in range(self.nPart)]
        for p in range(self.nPart):
            posList = self.truePositions[p]
            d = 0
            for s in range(1,len(posList)):
                d += np.linalg.norm(posList[s] - posList[s-1])
                if s in self.truePosStepIndices[p]:
                    stepLengths[p].append(d)
                    d = 0
        return stepLengths
