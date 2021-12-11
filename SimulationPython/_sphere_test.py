import numpy as np

from SimulationPython.Simulation.Sphere import Sphere

S = Sphere(0,0,0,-1,-1, 1, 1)

ray = np.array([np.array([0, 0, 0]), np.array([1, 0, 0])])
print(S.findIntersection(ray, 100))

ray = np.array([np.array([0, 0, 0]), np.array([-1, 0, 0])])
print(S.findIntersection(ray, 100))

ray = np.array([np.array([0, 0, 0]), np.array([1/2**0.5, 1/2**0.5, 0])])
print(S.findIntersection(ray, 100))

ray = np.array([np.array([0, 0, 0]), np.array([1, 0, 0])])
S.setPos(np.array([2.,0., 0.]))
print(S.findIntersection(ray, 100))

ray = np.array([np.array([0, 0, 0]), np.array([1, 0, 0])])
S.setPos(np.array([-2.,0., 0.]))
print(S.findIntersection(ray, 100))

ray = np.array([np.array([0, 0, 0]), np.array([-1, 0, 0])])
S.setPos(np.array([-2.,0., 0.]))
print(S.findIntersection(ray, 100))

ray = np.array([np.array([0, 4, 0]), np.array([1, 0, 0])])
print(S.findIntersection(ray, 100))
