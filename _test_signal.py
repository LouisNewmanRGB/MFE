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
'''
x = [0, 1, 2]
y = x
yerr = [1, 1, 1]
plt.errorbar(x, y, yerr, color="red", ecolor="green")
plt.grid()
plt.show()
'''


D = 2
diffusionTime = 2
radius = 8
print("start")
truePDF = Util.getPDF_sphere(D, diffusionTime, radius, 500, 5)
trueCDF = Util.getCDF_sphere(D, diffusionTime, radius, 500, 5)

pdfPointsX = np.linspace(0, 1.1*radius, 500)
#pdfPointsY = [truePDF(p) for p in pdfPointsX]
pdfPointsY = truePDF(pdfPointsX)
print("done")
