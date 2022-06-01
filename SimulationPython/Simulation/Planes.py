import numpy as np

from SimulationPython.Simulation.AbstractCompartment import AbstractCompartment

class Planes(AbstractCompartment):
    def __init__(self, T2, diffusivity, spacing):
        """Two parallel planes (perpendicular to x axis)"""
        super(Planes, self).__init__(None, None, None, T2, diffusivity)
        self.halfSpacing = spacing/2

    def findIntersection(self, ray, maxDistance):
        direction = ray[1]
        dir_fraction = np.inf
        if direction[0] != 0:
            dir_fraction = 1/direction[0]

        t1 = (self.halfSpacing - ray[0,0])*dir_fraction
        t2 = (-self.halfSpacing - ray[0,0])*dir_fraction

        if t2 > t1:
            t1, t2 = t2, t1
        if t1 > 0:
            return t1
        elif t2 > 0:
            return t2
        return np.inf

    def collide(self, particle, oldPos, intersection, sim):
        v = particle.getVelocity()
        v[0] = -v[0]
        particle.setVelocity(v)

    def contains(self, pos):
        return pos[0] > -self.halfSpacing and pos[0] < self.halfSpacing

    def setParentComp(self, comp):
        pass

    def plot(self, ax):
        pass
