import scipy.io
import numpy as np

from Util import Util

radius = 8
diffusionTimes = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 20, 25, 30, 40, 50, 60, 75, 100]
D = 2
qPoints = np.linspace(0, 0.4, 101)[1:]

#jnpZeroArray = Util.Jnp_zeros(100,100, 10)
#scipy.io.savemat("jnpZeroArray.mat", dict(jnpZeroArray=jnpZeroArray))


data = np.empty((len(diffusionTimes), len(qPoints)))
for t in range(len(diffusionTimes)):
    dt = diffusionTimes[t]
    signalFunc = Util.getSignal_sphere_fin(radius, D, dt, 100, 100, 10)
    data[t,:] = signalFunc(qPoints)

scipy.io.savemat("signals.mat", dict(data=data))
