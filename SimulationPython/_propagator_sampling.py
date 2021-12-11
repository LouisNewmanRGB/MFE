import numpy as np
import random
import matplotlib.pyplot as plt
import scipy.stats

from SimulationPython.Simulation.Util import Util

a = 10
D = 2
t = 4

def propa(z0, z1):
    return 1/a + (2/a)*np.sum([np.exp(-n**2 * np.pi**2 *D*t/a**2)*np.cos(n*np.pi*z0/a)*np.cos(n*np.pi*z1/a) for n in range(1, 11)], axis=0)

#points = np.linspace(0, a, 500)
#plt.plot(points, propa(0, points))
#plt.show()

zNumber = 100000

disp = np.empty(zNumber)

for i in range(zNumber):
    z0 = Util.getRandomU(a) + a/2 ###########
    next = False
    while not(next):
        z1 = Util.getRandomU(a) + a/2
        proba = propa(z0, z1)
        assert proba <= 1
        if random.random() < proba:
            disp[i] = z1 - z0
            next = True

bw = 2*scipy.stats.iqr(disp, rng=(25, 75))/(len(disp))**(1/3)
nBins = int((np.max(disp) - np.min(disp))/bw)
plt.hist(disp, bins=nBins, density = True, stacked=True)
plt.show()

def signal(q):
    return np.transpose(np.abs( np.average(np.exp(1j*np.matmul(np.transpose(np.matrix(disp)), np.matrix(q))), axis=0) ))

theoretical = Util.getSignal_plane_fin(a, D, t, 10)
qPoints = np.linspace(0, 1.5, 101)[1:]
plt.plot(qPoints, signal(qPoints))
plt.plot(qPoints, theoretical(qPoints), color="red")
plt.yscale("log")
plt.show()

