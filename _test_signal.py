import matplotlib.pyplot as plt
import numpy as np
import scipy.special

from Util import Util

radius = 8
nNumber = 4
nt = 5
pointsTest = np.linspace(0, 21, 500)[1:]

for n in range(nNumber):
    plt.plot(pointsTest, scipy.special.spherical_jn(n, pointsTest, derivative=True))
#plt.plot(pointsTest, scipy.special.spherical_jn(n, pointsTest, derivative=False), color="red")
#xPoints = (2.*np.arange(1,nt+n+1)-1)*np.pi/2
plt.legend(list(range(nNumber)))
xPoints2 = np.arange(1,nt+nNumber+1)*np.pi
xPoints = Util.genAlphaList(nt, nt, radius)*radius
plt.scatter(xPoints, np.zeros(len(xPoints)), color="black")
#plt.scatter(xPoints2, np.zeros(len(xPoints2)), color="red")
plt.grid()
plt.show()

points = np.linspace(0, 1.5, 500)[1:]
radius = 8
D = 2
t = 20
nNumber = 10
zeroNumber = 10

qStarPoints = np.linspace(0, 3.6, 500)[1:]
for tStar in [0.05, 0.1, 0.2, 0.5, 1]:
    fStar = Util.getSignal_sphere_red(tStar, nNumber, zeroNumber)
    #fStar = Util.getSignal_cylinder(tStar, nNumber, zeroNumber)
    #fStar = Util.getSignal_plane(tStar, 100)
    plt.plot(qStarPoints, fStar(qStarPoints))
plt.yscale("log")
plt.grid()
plt.show()


f = Util.getSignal_sphere_inf(radius)
g = Util.getSignal_sphere_fin(radius, D, t, nNumber, zeroNumber)
#h = Util.getSignal_sphere_2(radius, D, t, 10*nNumber, 10*zeroNumber)

plt.plot(points, f(points), color="red")
plt.plot(points, g(points), color="green")
#plt.plot(points, h(points), color="blue")
plt.yscale("log")
plt.show()
