import numpy as np
import pyrr.geometric_tests as gt #installed for ray calculations

from AbstractCompartment import AbstractCompartment
from Simulation import Simulation

class Environment(AbstractCompartment):
    def __init__(self, T2, diffusivity, sizeX, sizeY, sizeZ):
        super(Environment, self).__init__(None, None, None, T2, diffusivity)
        self.size = np.array([sizeX, sizeY, sizeZ])
        self.tolerance = self.size*Simulation.TOL
        self.aabb = np.array([-self.size/2, self.size/2])

    def findIntersection(self, ray, maxDistance):
        #if (ray[0] >= self.size/2 - maxDistance).any() or (ray[0] <= self.size/2 + maxDistance).any():
        return gt.ray_intersect_aabb(ray, self.aabb)

    def collide(self, particle, oldPos, intersection, sim):
        line = intersection + particle.getVelocity()/particle.getSpeed()*self.tolerance
        newPos = intersection.copy()
        for i in range(3):
            if line[i] > self.size[i]/2:
                newPos[i] = -self.size[i]/2
            elif line[i] < -self.size[i]/2:
                newPos[i] = self.size[i]/2
        particle.setPos(newPos)
        newComp = sim.findCompartment(particle)
        particle.changeCompartment(newComp, sim.getTimeStep())

    def contains(self, particle):
        return True

    def plot(self, ax):
        pass

    def getSize(self):
        return self.size
