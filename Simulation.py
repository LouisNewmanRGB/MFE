import numpy as np
import random
import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm #color map

class Simulation():
    TIME_TOL = 1e-10
    
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
        
        #for trajectory plotting
        self.positions = None
        self.truePositions = None
        
    def getTimeStep(self):
        return self.timeStep
        
    def getPositions(self):
        return self.positions
    
    def getTruePositions(self):
        return self.truePositions
        
    def getRandomDirection(self):
        theta = random.random() * 2*math.pi
        phi = random.random() * math.pi
        return np.array([math.sin(phi) * math.cos(theta), math.sin(phi) * math.sin(theta), math.cos(phi)])
        
    def findCompartment(self, particle):
        for compartment in self.compartments:
            if compartment.contains(particle):
                return compartment
                break
    
    def run(self, seed=None, calcData = False):
        #calcData determines whether intermediate positions will be stored in memory
        
        nPart = len(self.particles)
        if seed !=None:
            random.seed(seed)
        if calcData:
            self.positions = np.zeros((nPart,self.nStep,3))
            self.truePositions = np.zeros((nPart,self.nStep,3))
        
        for p in range(nPart):
            particle = self.particles[p]
            particle.setVelocity(self.getRandomDirection())
            newComp = self.findCompartment(particle)
            particle.changeCompartment(newComp, self.timeStep)
            for n in range(self.nStep):
                if calcData:
                    self.positions[p,n,:] = particle.getPos()
                    self.truePositions[p,n,:] = particle.getTruePos()
                self.nextStep(particle)
                
        if calcData:
            self.positions = self.positions.transpose(1, 0, 2)
            self.truePositions = self.truePositions.transpose(1, 0, 2)
            
        self.calcDisplacements()
        self.calcDistances()
    
    def nextStep(self, particle):
        t = 0 #elapsed time during the step
        particle.setVelocity(particle.getSpeed()*self.getRandomDirection()) #random direction at each step
        
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
                compartment.collide(particle, intersections[firstIndex], reachTime, self)
                t += reachTime
            else:
                particle.move(self.timeStep - t)
                t = self.timeStep
        
    def plot(self, positionType = True):
        #3D plot of trajectories
        #positionType represents whether to plot "true" or "in cell" positions
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        nPart = len(self.particles)
        colors = cm.rainbow(np.linspace(0, 1, nPart))
        
        if positionType:
            data = self.truePositions
        else:
            ax.set_xlim3d([self.environment.getPos()[0] - self.environment.getSize()[0]/2,
                            self.environment.getPos()[0] + self.environment.getSize()[0]/2])
            ax.set_ylim3d([self.environment.getPos()[1] - self.environment.getSize()[1]/2,
                            self.environment.getPos()[1] + self.environment.getSize()[1]/2])
            ax.set_zlim3d([self.environment.getPos()[2] - self.environment.getSize()[2]/2,
                            self.environment.getPos()[2] + self.environment.getSize()[2]/2])
            data = self.positions
        
        for p in range(nPart):
            xArray = data[:, p, 0]
            yArray = data[:, p, 1]
            zArray = data[:, p, 2]
            ax.scatter(xArray, yArray, zArray, color = colors[p])
            ax.plot(xArray, yArray, zArray, color = colors[p])
        
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        plt.show()
        
    def calcDisplacements(self):
        endPos = np.array([particle.getTruePos() for particle in self.particles])
        self.displacements = endPos - self.startingPositions
        
    def calcDistances(self):
        self.distances = np.array([np.linalg.norm(disp) for disp in self.displacements])
        
    def getDistances(self):
        return self.distances