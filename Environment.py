import numpy as np
import pyrr.geometric_tests as gt #installed for ray calculations

from AbstractCompartment import AbstractCompartment
from Simulation import Simulation

class Environment(AbstractCompartment):
    def __init__(self, x, y, z, T2, diffusivity, sizeX, sizeY, sizeZ):
        super(Environment, self).__init__(x, y, z, T2, diffusivity)
        self.size = np.array([sizeX, sizeY, sizeZ])

    def findIntersection(self, ray, maxDistance):
        boundingBox = np.array([self.pos - self.size/2, self.pos + self.size/2])
        return gt.ray_intersect_aabb(ray, boundingBox)

    def collide(self, particle, oldPos, intersection, sim):
        line = intersection + particle.getVelocity()*Simulation.TIME_TOL
        newPos = intersection.copy()
        for i in range(3):
            if line[i] > self.pos[i] + self.size[i]/2:
                newPos[i] = self.pos[i] - self.size[i]/2
            elif line[i] < self.pos[i] - self.size[i]/2:
                newPos[i] = self.pos[i] + self.size[i]/2
        particle.setPos(newPos)
        newComp = sim.findCompartment(particle)
        particle.changeCompartment(newComp, sim.getTimeStep())

    def contains(self, particle):
        return True

    def plot(self, ax):
        pass

    def getSize(self):
        return self.size
