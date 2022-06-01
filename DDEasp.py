import numpy as np
import matplotlib.pyplot as plt

class DDEasp():
    def __init__(self, delta, Delta, psi):
        self.delta = delta
        self.Delta = Delta
        self.psi = psi

    def getTimes(self):
        return np.array([self.delta, self.Delta, self.Delta + self.delta,
                         self.Delta + 2*self.delta, 2*self.Delta + self.delta, 2*self.Delta + 2*self.delta,
                         2*self.Delta + 3*self.delta])

    def getDirections(self):
        psi = self.psi*np.pi/180
        return np.array([[1, 0, 0], [0, 0, 0], [-1, 0, 0],
                         [-np.cos(psi), -np.sin(psi), 0], [0, 0, 0], [np.cos(psi), np.sin(psi), 0],
                         [0, 0, 0]])

    def plot(self):
        xMaxs = self.getTimes()
        xMins = [0] + list(xMaxs[:-1])
        ys = [dir[0] for dir in self.getDirections()]
        plt.hlines(ys, xmin=xMins, xmax=xMaxs)
        plt.scatter(xMaxs, ys)
        plt.xlabel("time [ms]")
        plt.ylabel("x Component of magnetic field gradient [au]")
        plt.title("DDE sequence of parameters \n"
                  "Delta={Delta}ms, delta={delta}ms, psi={psi}deg".format(Delta=self.Delta, delta=self.delta, psi=self.psi))
        plt.grid()
        plt.show()
