import numpy as np

from AbstractCompartment import AbstractCompartment
from Simulation import Simulation

class Environment(AbstractCompartment):
    def __init__(self, T2, diffusivity, sizeX, sizeY, sizeZ):
        super(Environment, self).__init__(None, None, None, T2, diffusivity)
        self.size = np.array([sizeX, sizeY, sizeZ])
        self.tolerance = max(self.size)*Simulation.TOL
        self.aabb = np.array([-self.size/2, self.size/2])

    def findIntersection(self, ray, maxDistance):
        #if (ray[0] >= self.size/2 - maxDistance).any() or (ray[0] <= self.size/2 + maxDistance).any():
        #return gt.ray_intersect_aabb([ray.origin, ray.direction], self.aabb)

        direction = ray[1]
        dir_fraction = np.empty(3, dtype = ray.dtype)
        dir_fraction[direction == 0.0] = np.inf
        dir_fraction[direction != 0.0] = np.divide(1.0, direction[direction != 0.0])

        t1 = (self.aabb[0,0] - ray[0,0]) * dir_fraction[ 0 ]
        t2 = (self.aabb[1,0] - ray[0,0]) * dir_fraction[ 0 ]
        t3 = (self.aabb[0,1] - ray[0,1]) * dir_fraction[ 1 ]
        t4 = (self.aabb[1,1] - ray[0,1]) * dir_fraction[ 1 ]
        t5 = (self.aabb[0,2] - ray[0,2]) * dir_fraction[ 2 ]
        t6 = (self.aabb[1,2] - ray[0,2]) * dir_fraction[ 2 ]

        tmin = max(min(t1, t2), min(t3, t4), min(t5, t6))
        tmax = min(max(t1, t2), max(t3, t4), max(t5, t6))

        if tmax < 0:
            return np.inf
        if tmin > tmax:
            return np.inf
        t = min(x for x in [tmin, tmax] if x >= 0)
        return t

    def collide(self, particle, oldPos, intersection, sim):
        line = intersection + particle.getVelocity()/particle.getSpeed()*self.tolerance
        newPos = intersection.copy()
        for i in range(3):
            if line[i] > self.size[i]/2:
                newPos[i] = -self.size[i]/2
            elif line[i] < -self.size[i]/2:
                newPos[i] = self.size[i]/2
        particle.setPos(newPos)
        newComp = sim.findCompartment(particle.getPos())
        particle.changeCompartment(newComp, sim.getTimeStep())

    def contains(self, particle):
        return True

    def plot(self, ax):
        pass

    def getSize(self):
        return self.size
