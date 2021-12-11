import matplotlib.pyplot as plt
import numpy as np

from SimulationPython.Simulation.Util import Util
"""
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
"""
print(help(Util.getSignal_cylinder_red))

#White et Dale
D = 1
t = 60
radius = 5
bValue = 4
print("White et Dale maximum qR value:", radius*(bValue/t)**0.5)

#GPU article
radius = 2
D = 0.5
tList = [50, 250]
nNumber = 10
zeroNumber = 10
bValue = 60

bPoints = np.linspace(0, bValue, 501)[1:]
for t in tList:
    signalFunc = Util.getSignal_cylinder_fin(radius, D, t, 10, 10)
    plt.plot(bPoints, signalFunc((bPoints/t)**0.5))
    print("GPU maximum qR value:", radius*(bPoints[-1]/t)**0.5)
print("My graphs go to qR=12")
plt.yscale("log")
plt.show()

"""
qStarPoints = np.linspace(0, 3.6, 500)[1:]
for tStar in [0.05, 0.1, 0.2, 0.5, 1]:
    #fStar = Util.getSignal_sphere_red(tStar, nNumber, zeroNumber)
    fStar = Util.getSignal_cylinder(tStar, nNumber, zeroNumber)
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
"""
