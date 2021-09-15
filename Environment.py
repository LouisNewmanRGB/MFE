import numpy as np
import pyrr.geometric_tests as gt #installed for ray calculations

from AbstractCompartment import AbstractCompartment
from Simulation import Simulation

class Environment(AbstractCompartment):
    def __init__(self, x, y, z, T2, diffusivity, sizeX, sizeY, sizeZ):
        super(Environment, self).__init__(x, y, z, T2, diffusivity)
        self.size = np.array([sizeX, sizeY, sizeZ])

    def findIntersection(self, particle):
        ray = np.array([particle.getPos() + particle.getVelocity()*Simulation.TIME_TOL, particle.getVelocity()])
        boundingBox = np.array([self.pos - self.size/2, self.pos + self.size/2])
        return gt.ray_intersect_aabb(ray, boundingBox)

    def collide(self, particle, intersection, reachTime, sim):
        newPos = 2 * intersection - particle.getPos() #position outside the bounding box
        for i in range(3):
            if newPos[i] > self.pos[i] + self.size[i]/2:
                newPos[i] = self.pos[i] - self.size[i]/2
            elif newPos[i] < self.pos[i] - self.size[i]/2:
                newPos[i] = self.pos[i] + self.size[i]/2
        particle.move(reachTime)
        particle.setPos(newPos)
        newComp = sim.findCompartment(particle)
        particle.changeCompartment(newComp, sim.getTimeStep())

    def contains(self, particle):
        return True

    def getSize(self):
        return self.size