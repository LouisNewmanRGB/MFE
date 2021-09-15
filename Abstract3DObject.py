import abc #for abstract classes
import numpy as np

class Abstract3DObject(metaclass=abc.ABCMeta):
   @abc.abstractmethod
   def __init__(self, x, y, z):
       self.pos = np.array([z, y, z]) #position of the 3D object in space

   def getPos(self):
       return self.pos

   def setPos(self, pos):
       self.pos = pos