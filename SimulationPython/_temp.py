import numpy as np
import matplotlib.pyplot as plt

from SimulationPython.Simulation.Util import Util
from SimulationPython.Simulation.StejskalTanner import StejskalTanner
from SimulationPython.Simulation.Particle3DPiecewiseGradient import Particle3DPiecewiseGradient
from SimulationPython.Simulation.Environment import Environment

t = 10
TE = 70
delta = 20
Delta = 30
G = 0.1
#W&D
TE = 140
delta = 60
Delta = delta
t = TE/2 - Delta
seq = StejskalTanner(t, TE, delta, Delta)

nStep = 3
timeSteps = [TE/nStep]*nStep
T2 = np.inf
env = Environment(T2, 2, 5, 5, 5)
part = Particle3DPiecewiseGradient(0, 0, 0, seq)
part.setVelocity(np.array([1, 0, 0]))
part.changeCompartment(env, TE/nStep)
for ts in timeSteps:
    part.move(ts)

seq.plot()
stop

"""
import numpy as np
import matplotlib.pyplot as plt

from Util import Util

points = np.linspace(0, 50, 50)
D = 2
a = 8

data = np.empty(len(points))
for t in range(len(points)):
    data[t] = Util.RMSD_planes(a, D, points[t], 20)

plt.plot(points, data)
plt.grid()
plt.show()
"""

a = 8
#alphaList = np.array([0, 4.4934, 7.7253, 10.9041, 14.0662, 17.2208])/a

alphaList = [(2*k + 1)*np.pi/2 for k in range(500)]
for j in range(5):
    for k in range(len(alphaList)):
        alphaList[k] = np.arctan(alphaList[k]) + k*np.pi
alphaList = alphaList[1:]
print(alphaList)

#alphaList = np.array([4.4934, 7.7253, 10.9041, 14.0662, 17.2208])/a
D = 2
t = 10

def f(r):
    #return 3/(4*np.pi*a**3) + 2*np.sum([np.exp(-D*t*alpha**2)*np.sin(alpha*r)*alpha/(4*np.pi*np.sin(alpha*a)**2) for alpha in alphaList])/(a*r)
    return 3*r**2/(a**3) + 2*r*np.sum([np.exp(-D*t*alpha**2)*np.sin(alpha*r)*alpha/(np.sin(alpha*a)**2) for alpha in alphaList])/a

def h(r):
    return 3*r**2/(a**3)

g = Util.getPDF_sphere(D, t, a, 1000, 5)
G = Util.getCDF_sphere(D, t, a, 500, 5)

points = np.linspace(0, 1.1*a, 500)
plt.plot(points, [g(r) for r in points])
plt.plot(points, [h(r) for r in points])
plt.grid()
plt.show()

plt.plot(points, G(points))
plt.plot(points, [G(r) for r in points])
plt.grid()
plt.show()
