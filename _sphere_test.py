import numpy as np

from Sphere import Sphere
from Particle3D import Particle3D
from Util import Util

S = Sphere(0,0,0,-1,-1, 1)
P = Particle3D(0.,0.,0.)

P.setVelocity(np.array([0.5,0., 0.]))
print(S.findIntersection(P))

P.setVelocity(np.array([-0.5,0., 0.]))
print(S.findIntersection(P))

P.setVelocity(np.array([0.5,0.5, 0.]))
print(S.findIntersection(P))

P.setVelocity(np.array([0.5,0., 0.]))
S.setPos(np.array([2.,0., 0.]))
print(S.findIntersection(P))

P.setVelocity(np.array([0.5,0., 0.]))
S.setPos(np.array([-2.,0., 0.]))
print(S.findIntersection(P))

P.setVelocity(np.array([-0.5,0., 0.]))
S.setPos(np.array([-2.,0., 0.]))
print(S.findIntersection(P))

P.setPos(np.array([0., 4., 0.]))
print(S.findIntersection(P))

print(Util.recursiveMax([[0,1],[2,-3]]))
