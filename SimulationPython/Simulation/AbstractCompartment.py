import abc #for abstract classes

from SimulationPython.Simulation.Abstract3DObject import Abstract3DObject

class AbstractCompartment(Abstract3DObject):
    def __init__(self, x, y, z, T2, diffusivity):
        super(AbstractCompartment, self).__init__(x, y, z)
        self.T2 = T2
        self.diffusivity = diffusivity

    @abc.abstractmethod
    def contains(self, particle):
       pass

    @abc.abstractmethod
    def findIntersection(self, ray, maxDistance):
       #cast a ray and find a possible intersection between ray and compartment
       pass

    @abc.abstractmethod
    def collide(self, particle, oldPos, intersection, sim):
       #collide with a particle at position intersection and time reachTime in simulation sim
       pass

    @abc.abstractmethod
    def plot(self, ax):
       #plot compartment on matplotlib 3D axes ax
       pass

    def getT2(self):
        return self.T2

    def getDiffusivity(self):
        return self.diffusivity
