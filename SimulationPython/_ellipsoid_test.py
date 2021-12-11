import numpy as np

from SimulationPython.Simulation.Ellipsoid import Ellipsoid
from SimulationPython.Simulation.Sphere import Sphere
from SimulationPython.Simulation.Particle3D import Particle3D

#E = Ellipsoid(0, 0, 0, -1, -1, -1, np.array([1, 1, 0]), np.array([0, 2, 2]), np.array([0.5, 0, 0.5]))
E = Ellipsoid(0, 0, 0, -1, -1, -1, 1, 1, 1)
S = Sphere(0,0,0,-1,-1, -1, 1)

ray = np.array([[0, 0, 0], [1, 0, 0]])
print(S.findIntersection(ray, 10))
print(E.findIntersection(ray, 10), "\n")

ray = np.array([[0, 0, 0], [-1, 0, 0]])
print(S.findIntersection(ray, 10))
print(E.findIntersection(ray, 10), "\n")

ray = np.array([[0, 0, 0], [2**-0.5, 2**-0.5, 0]])
print(S.findIntersection(ray, 10))
print(E.findIntersection(ray, 10), "\n")

ray = np.array([[0, 0, 0], [1, 0, 0]])
S.setPos(np.array([2.,0., 0.]))
E.setPos(np.array([2.,0., 0.]))
print(S.findIntersection(ray, 10))
print(E.findIntersection(ray, 10), "\n")

ray = np.array([[0, 0, 0], [1, 0, 0]])
S.setPos(np.array([-2.,0., 0.]))
E.setPos(np.array([-2.,0., 0.]))
print(S.findIntersection(ray, 10))
print(E.findIntersection(ray, 10), "\n")

ray = np.array([[0, 0, 0], [-1, 0, 0]])
S.setPos(np.array([-2.,0., 0.]))
E.setPos(np.array([-2.,0., 0.]))
print(S.findIntersection(ray, 10))
print(E.findIntersection(ray, 10), "\n")

ray = np.array([[0, 4, 0], [-1, 0, 0]])
print(S.findIntersection(ray, 10))
print(E.findIntersection(ray, 10), "\n")



P = Particle3D(0, 0, 0)
P.setPos(np.array([0, 0, 0]))
E = Ellipsoid(0, 0, 0, -1, -1, -1, 0.5, 2, 1)
print(E.contains(P))
P.setPos(np.array([0.6, 0, 0]))
print(E.contains(P))
P.setPos(np.array([0, 2.5, 0]))
print(E.contains(P))
P.setPos(np.array([0, 0, .99]))
print(E.contains(P))
