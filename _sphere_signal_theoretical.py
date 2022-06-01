import matplotlib.pyplot as plt
import numpy as np

from Util import Util

D = 2
r = 8
DeltaList = [10, 20, 35, 50]

q = np.linspace(0, 1.15, 501)[1:]

for Delta in DeltaList:
    s = Util.getSignal_sphere_fin(r, D, Delta, 20, 20, 5)
    plt.plot(q, s(q))

sInf = Util.getSignal_sphere_inf(r)
plt.plot(q, sInf(q))
plt.legend(["Diffusion time ={Delta}ms".format(Delta=Delta) for Delta in DeltaList] + ["Infinite diffusion time"])
plt.yscale("log")
plt.ylim(1e-7)
plt.xlabel("q=gamma*G*delta [um-1]")
plt.ylabel("Signal attenuation")
plt.title("Theoretical SGP Signal Attenuation for Diffusion in an Impermeable Sphere")
plt.grid()
plt.show()
