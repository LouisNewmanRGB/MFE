import numpy as np

from Sphere import Sphere
from Util import Util

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

#print(Util.recursiveMax([[0,1],[2,-3]]))
