import numpy as np
import random

from SimulationPython.Simulation.AbstractCompartment import AbstractCompartment


class Ellipsoid(AbstractCompartment):
    def __init__(self, x, y, z, T2, diffusivity, permeability, a, b, c):
        super(Ellipsoid, self).__init__(x, y, z, T2, diffusivity)
        self.permeability = permeability
        self.parentComp = None
        if type(a) != np.ndarray:
            self.R = np.diag([a, b, c])
        else:
            self.R = np.zeros((3,3))
            self.R[:,0], self.R[:,1], self.R[:,2] = a, b, c
        self.maxRadius = np.max([np.linalg.norm(a), np.linalg.norm(b), np.linalg.norm(c)])
        self.invR = np.linalg.inv(self.R)
        self.A = np.matmul(np.transpose(self.invR), self.invR) #parameter in quadratic form

    def setParentComp(self, compartment):
        self.parentComp = compartment

    def findIntersection(self, ray, maxDistance):
        if np.linalg.norm(ray[0] - self.pos) - self.maxRadius <= maxDistance:
            L = np.matmul(self.invR,(self.pos - ray[0]))
            invRD = np.matmul(self.invR,ray[1])
            a = np.linalg.norm(invRD)**2
            b = -2*np.dot(invRD, L)
            c = np.linalg.norm(L)**2 - 1
            delta = b**2 - 4*a*c
            if delta > 0:
                t = (-b - delta**0.5)/(2*a)
                if t > 0:
                    return t
                else:
                    t = (-b + delta**0.5)/(2*a)
                    if t > 0:
                        return t
            elif delta == 0:
                t = -b/(2*a)
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
            normal = 2*np.matmul(self.A, intersection - self.pos)
            normal = normal/np.linalg.norm(normal)
            particle.setVelocity(particle.getVelocity() - 2*normal*np.dot(normal, particle.getVelocity()))
        else:
            particle.changeCompartment(otherSideComp, sim.getTimeStep())

    def contains(self, pos):
        return np.linalg.norm(np.matmul(self.invR, pos - self.pos)) < 1

    def plot(self, ax):
        ############################################
        u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
        t = np.array([np.cos(u)*np.sin(v), np.sin(u)*np.sin(v), np.cos(v)])
        p = np.zeros(t.shape)
        for i in range(t.shape[1]):
            for j in range(t.shape[2]):
                p[:,i,j] = self.pos + np.matmul(self.R, t[:,i,j])
        ax.plot_wireframe(p[0,:,:], p[1,:,:], p[2,:,:], color="black")
        #############################################
