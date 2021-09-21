import numpy as np
import random

from AbstractCompartment import AbstractCompartment
from Simulation import Simulation
from Util import Util
from Particle3D import Particle3D

class Sphere(AbstractCompartment):
    def __init__(self, x, y, z, T2, diffusivity, probInOut, radius):
        super(Sphere, self).__init__(x, y, z, T2, diffusivity)
        self.probInOut = probInOut
        self.radius = radius

    def findIntersection(self, particle):
        ray = np.array([particle.getPos() + particle.getVelocity()*Simulation.TIME_TOL, particle.getVelocity()/particle.getSpeed()])
        L = self.pos - ray[0]
        tList = Util.rootsReal([1, -2*np.dot(ray[1], L), np.linalg.norm(L)**2 - self.radius**2])
        if len(tList) > 0:
            tList = np.sort(tList)
            if tList[0] > 0:
                return ray[0] + tList[0] * ray[1]
            elif len(tList) == 2 and tList[1] > 0:
                return ray[0] + tList[1] * ray[1]
        return None

    def collide(self, particle, oldPos, intersection, sim):
        otherSideComp = sim.findCompartment(Particle3D(*intersection + particle.getVelocity()*Simulation.TIME_TOL))
        if otherSideComp.contains(Particle3D(*oldPos)):
            #the particle is leaving a compartment
            probability = self.probInOut
        else:
            #the particle is entering a compartment
            probability = self.probInOut * (self.diffusivity/particle.getCompartment().getDiffusivity())**0.5
        if random.random() > probability:
            #deflection
            normal = (intersection - self.pos)
            normal = normal/np.linalg.norm(normal)
            particle.setVelocity(particle.getVelocity() - 2*normal*np.dot(normal, particle.getVelocity()))
        else:
            particle.changeCompartment(otherSideComp, sim.getTimeStep())

    def contains(self, particle):
        return np.linalg.norm(particle.getPos() - self.pos) < self.radius

    def plot(self, ax):
        ############################################
        u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
        x = self.pos[0] + self.radius*np.cos(u)*np.sin(v)
        y = self.pos[1] + self.radius*np.sin(u)*np.sin(v)
        z = self.pos[2] + self.radius*np.cos(v)
        ax.plot_wireframe(x, y, z, color="black")
        #############################################
