import numpy as np
import math

from Abstract3DObject import Abstract3DObject

class Particle3D(Abstract3DObject):
   def __init__(self, x, y, z, vx, vy, vz):
       super(Particle3D, self).__init__(float(x), float(y), float(z))
       self.truePos = self.pos.copy() #position without wrap-around effect
       self.velocity = np.array([vx, vy, vz])
       self.compartment = None
       self.signal = 1

   def __init__(self, x, y, z):
       super(Particle3D, self).__init__(float(x), float(y), float(z))
       self.truePos = self.pos.copy() #position without wrap-around effect
       self.velocity = None
       self.compartment = None
       self.signal = 1

   def move(self, time):
       self.pos += time*self.velocity
       self.truePos += time*self.velocity
       self.signal *= math.exp(-time/self.compartment.getT2()) #T2 signal attenuation

   def changeCompartment(self, compartment, timeStep):
       self.compartment = compartment
       self.changeSpeed((6*self.compartment.getDiffusivity()/timeStep)**0.5)

   def getCompartment(self):
        return self.compartment

   def getTruePos(self):
        return self.truePos

   def setVelocity(self, velocity):
       self.velocity = velocity

   def getVelocity(self):
       return self.velocity

   def getSpeed(self):
       return np.linalg.norm(self.velocity)

   def changeSpeed(self, newSpeed):
       self.velocity = newSpeed*self.velocity/self.getSpeed()

   def getSignal(self):
       return self.signal
