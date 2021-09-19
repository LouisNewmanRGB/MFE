import numpy as np
from statsmodels.distributions.empirical_distribution import ECDF
import matplotlib.cm as cm
import matplotlib.pyplot as plt
"""
data = np.random.normal(0,5, size=20)

ecdf = ECDF(data)
print(len(ecdf.x))
print(len(ecdf.y))
print(ecdf.x)
print(ecdf.y)
for i in range(len(ecdf.x)-1):
    plt.hlines(ecdf.y[i], ecdf.x[i], ecdf.x[i+1])

plt.show()
#plt.plot(ecdf.x,ecdf.y)
"""
"""
    #infinite norm error
    ecdf = ECDF(sim.getDistances())
    #print(ecdf.x)
    #print(ecdf.y)
    print(np.max(np.abs(ecdf.y[1:] - trueDistrib(ecdf.x[1:]))))
"""
diffusionTimes = [2, 4]
nStep = [1, 2]
errors = np.array([[0, 1], [2, 3]])
totalSteps=1e5

colors = cm.rainbow(np.linspace(0, 1, len(diffusionTimes)))
for t in range(len(diffusionTimes)):
    plt.scatter(nStep, errors[t,:], color = colors[t])
plt.xscale("log")
plt.xlabel("Number of steps")
plt.ylabel("Supremum of distances between exact and empirical cumulative distribution functions")
plt.legend(["Diffusion time = {diffusionTime}ms".format(diffusionTime=dt) for dt in diffusionTimes])
plt.title("Random walk simulation of free diffusion for different diffusion times and numbers of steps\n"
          "particles x steps = {totalSteps}".format(totalSteps=totalSteps))
plt.show()
