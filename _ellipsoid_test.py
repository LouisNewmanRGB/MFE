import numpy as np

from Ellipsoid import Ellipsoid
from Sphere import Sphere
from Particle3D import Particle3D
from Util import Util

E = Ellipsoid(0, 0, 0, -1, -1, -1, np.array([1, 1, 0]), np.array([0, 2, 2]), np.array([0.5, 0, 0.5]))
S = Sphere(0,0,0,-1,-1, -1, 0.5)
P = Particle3D(0.,0.,0.)

P.setVelocity(np.array([0.5,0., 0.]))
print(S.findIntersection(P))
print(E.findIntersection(P), "\n")

P.setVelocity(np.array([-0.5,0., 0.]))
print(S.findIntersection(P))
print(E.findIntersection(P), "\n")

P.setVelocity(np.array([0.5,0.5, 0.]))
print(S.findIntersection(P))
print(E.findIntersection(P), "\n")

P.setVelocity(np.array([0.5,0., 0.]))
S.setPos(np.array([2.,0., 0.]))
E.setPos(np.array([2.,0., 0.]))
print(S.findIntersection(P))
print(E.findIntersection(P), "\n")

P.setVelocity(np.array([0.5,0., 0.]))
S.setPos(np.array([-2.,0., 0.]))
E.setPos(np.array([-2.,0., 0.]))
print(S.findIntersection(P))
print(E.findIntersection(P), "\n")

P.setVelocity(np.array([-0.5,0., 0.]))
S.setPos(np.array([-2.,0., 0.]))
E.setPos(np.array([-2.,0., 0.]))
print(S.findIntersection(P))
print(E.findIntersection(P), "\n")

P.setPos(np.array([0., 4., 0.]))
print(S.findIntersection(P))
print(E.findIntersection(P), "\n")



#print(Util.recursiveMax([[0,1],[2,-3]]))
P.setPos(np.array([0, 0, 0]))
E = Ellipsoid(0, 0, 0, -1, -1, -1, 0.5, 2, 1)
print(E.contains(P))
P.setPos(np.array([0.6, 0, 0]))
print(E.contains(P))
P.setPos(np.array([0, 2.5, 0]))
print(E.contains(P))
P.setPos(np.array([0, 0, .99]))
print(E.contains(P))

E.plot(1)
