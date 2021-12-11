import numpy as np
import matplotlib.pyplot as plt

from SimulationPython.Simulation.Util import Util

nDir = 2000
directions = np.empty((nDir, 3))
for n in range(nDir):
    directions[n] = Util.getRandomDirection()

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
#ax.scatter(directions[:, 0], directions[:, 1], directions[:, 2]*5/4)
ax.scatter(directions[:, 1], directions[:, 2], directions[:, 0]*5/4)
ax.set_xlim3d([-1.1, 1.1])
ax.set_ylim3d([-1.1, 1.1])
ax.set_zlim3d([-1.1, 1.1])
plt.show()
