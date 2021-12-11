import numpy as np
import random

from SimulationPython.Simulation.AbstractCompartment import AbstractCompartment

class Sphere(AbstractCompartment):
    def __init__(self, x, y, z, T2, diffusivity, permeability, radius):
        super(Sphere, self).__init__(x, y, z, T2, diffusivity)
        self.permeability = permeability
        self.radius = radius
        self.parentComp = None

    def setParentComp(self, compartment):
        self.parentComp = compartment

    def findIntersection(self, ray, maxDistance):
        #ray direction vector must be normalized!
        if np.linalg.norm(ray[0] - self.pos) - self.radius <= maxDistance:
            L = self.pos - ray[0]
            b = -2*np.dot(ray[1], L)
            c = np.linalg.norm(L)**2 - self.radius**2
            delta = b**2 - 4*c
            if delta > 0:
                t = (-b - delta**0.5) / 2
                if t > 0:
                    return t
                else:
                    t = (-b + delta**0.5) / 2
                    if t > 0:
                        return t
            elif delta == 0:
                t = -b/2
                if t > 0:
                    return t
        return np.inf

    def collide(self, particle, oldPos, intersection, sim):
        if self.contains(oldPos):
            #we are leaving this compartment
            otherSideComp = self.parentComp
        else:
            #we are entering this compartment
            otherSideComp = self
        transmissionProba = 4*self.permeability/particle.getSpeed()
        assert(transmissionProba <= 1)

        if random.random() > transmissionProba:
            #deflection
            normal = (intersection - self.pos)
            normal = normal/np.linalg.norm(normal)
            particle.setVelocity(particle.getVelocity() - 2*normal*np.dot(normal, particle.getVelocity()))
        else:
            particle.changeCompartment(otherSideComp, sim.getTimeStep())

    def contains(self, pos):
        return np.linalg.norm(pos - self.pos) < self.radius

    def plot(self, ax):
        ############################################
        u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
        x = self.pos[0] + self.radius*np.cos(u)*np.sin(v)
        y = self.pos[1] + self.radius*np.sin(u)*np.sin(v)
        z = self.pos[2] + self.radius*np.cos(v)
        ax.plot_wireframe(x, y, z, color="black")
        #############################################
