import matplotlib.pyplot as plt
import numpy as np

from Util import Util

"""
points = np.linspace(0, 1.5, 500)[1:]

#def f(qR):
#    return 9*(qR*np.cos(qR) - np.sin(qR))**2 / qR**6

f = Util.getSignal_sphere(8)

plt.plot(points, f(points))
plt.yscale("log")
plt.show()
"""

x = [0, 1, 2]
y = x
yerr = [1, 1, 1]
plt.errorbar(x, y, yerr, color="red", ecolor="green")
plt.grid()
plt.show()

