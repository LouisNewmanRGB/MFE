import numpy as np
import random

from AbstractCompartment import AbstractCompartment
from Simulation import Simulation
from Util import Util
from Particle3D import Particle3D

class Ellipsoid(AbstractCompartment):
    def __init__(self, x, y, z, T2, diffusivity, probInOut, a, b, c):
        super(Ellipsoid, self).__init__(x, y, z, T2, diffusivity)
        self.probInOut = probInOut
        if type(a) != np.ndarray:
            self.R = np.diag([a, b, c])
        else:
            self.R = np.zeros((3,3))
            self.R[:,0], self.R[:,1], self.R[:,2] = a, b, c
        self.maxRadius = np.max([np.linalg.norm(a), np.linalg.norm(b), np.linalg.norm(c)])
        self.invR = np.linalg.inv(self.R)
        self.A = np.matmul(np.transpose(self.invR), self.invR) #parameter in quadratic form

    def findIntersection(self, ray, maxDistance):
        if np.linalg.norm(ray[0] - self.pos) - self.maxRadius <= maxDistance:
            L = np.matmul(self.invR,(self.pos - ray[0]))
            invRD = np.matmul(self.invR,ray[1])
            tList = Util.rootsReal([np.linalg.norm(invRD)**2, -2*np.dot(invRD, L), np.linalg.norm(L)**2 - 1])
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
            normal = 2*np.matmul(self.A, intersection - self.pos)
            normal = normal/np.linalg.norm(normal)
            particle.setVelocity(particle.getVelocity() - 2*normal*np.dot(normal, particle.getVelocity()))
        else:
            particle.changeCompartment(otherSideComp, sim.getTimeStep())

    def contains(self, particle):
        return np.linalg.norm(np.matmul(self.invR, particle.getPos() - self.pos)) < 1

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
