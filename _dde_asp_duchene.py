import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
import scipy.integrate

from SimulationCpp import Simulation
from DDEasp import DDEasp
from SimulationCppResults import SimulationCppResults

simulateDict = False
nStep = 50
nPart = int(5*1e5)
delta = 12 #ms
Delta = 60 #ms
diffusionTime = 2*Delta + 2*delta
T2 = np.inf
D = 2 #um2/ms

qValues = 1e-3 * np.array([7, 14, 21, 28, 35]) #mm-1
angles = [0, 180]

#data
experimentalSignal = np.array([1, 0.8, 0.84, 0.4, 0.51, 0.144, 0.248, 0.048, 0.112, 0.02, 0.056]) #values read on GD's graph
#Histology data
histologyRadiuses = [3, 6.5, 10, 13, 21, 24, 28, 35, 38, 41, 46, 49, 59]
histologyPSD = np.array([7, 26.5, 25.6, 14.2, 5.8, 4, 2.8, 0.8, 0.4, 0.3, 0.3, 0.2, 0])
histologyPSD /= np.trapz(histologyPSD, histologyRadiuses) #normalization

#dictionnary of 20 + 5 values to estimate
radiuses = np.linspace(3, 70, 20)
diffusionCoefs = np.linspace(1, 2, 5)
if simulateDict:
    dict = np.empty((len(angles)*len(qValues) + 1, len(radiuses) + len(diffusionCoefs)))
    dictStds = np.empty((len(angles)*len(qValues) + 1, len(radiuses) + len(diffusionCoefs)))
    for p in range(len(angles)):
        psi = angles[p]
        for r in range(len(radiuses)):
            radius = radiuses[r]
            envSize = 5*radius
            nPartEff = nPart
            nStepEff = nStep
            sim = Simulation(nStepEff, diffusionTime/nStepEff, T2, D, envSize, envSize, envSize)
            sim.addCylinderBasic(0, 0, T2, D, 0, radius)
            sim.createStartingPositions(nPartEff, True)
            sequence = DDEasp(delta, Delta, psi)
            sim.createSequencePWG(sequence.getTimes(), sequence.getDirections())
            sim.run()
            dataVector = sim.getResults()
            results = SimulationCppResults(dataVector, nPartEff)

            dict[0, r], dictStds[0, r] = 1, 0
            for q in range(len(qValues)):
                dict[1+len(angles)*q + p, r], dictStds[1+len(angles)*q + p, r] = \
                    results.getFiniteGradientSignalGammaG(2*np.pi*qValues[q]/delta, includeStd=True)[0:2]

        for d in range(len(diffusionCoefs)):
            coef = diffusionCoefs[d]
            envSize = 10
            sim = Simulation(nStep, diffusionTime/nStep, T2, coef, envSize, envSize, envSize)
            sim.createStartingPositions(nPart, False)
            sequence = DDEasp(delta, Delta, psi)
            sim.createSequencePWG(sequence.getTimes(), sequence.getDirections())
            sim.run()
            dataVector = sim.getResults()
            results = SimulationCppResults(dataVector, nPart)

            dict[0, len(radiuses) + d], dictStds[0, len(radiuses) + d] = 1, 0
            for q in range(len(qValues)):
                dict[1+len(angles)*q + p, len(radiuses) + d], dictStds[1+len(angles)*q + p, len(radiuses) + d] = \
                    results.getFiniteGradientSignalGammaG(2*np.pi*qValues[q]/delta, includeStd=True)[0:2]

    np.save("numpy saves/AMICO_dict_duchene_500k.npy", dict)
else:
    dict = np.load("numpy saves/AMICO_dict_duchene_500k.npy")

#optimization: solving NNLS problem
n = len(radiuses) + len(diffusionCoefs)

#L-curve method -> find the "corner"
lambdas = [0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.1, 0.2, 0.5, 0.7, 0.9, 1, 2]
residualNorms = []
solutionNorms = []
for l in range(len(lambdas)):
    lam = lambdas[l]
    A = np.concatenate((dict, lam*np.identity(n)))
    b = np.concatenate((experimentalSignal, np.zeros(n)))
    estimation, r = scipy.optimize.nnls(A, b)

    residualNorms.append(np.linalg.norm(np.matmul(dict, estimation) - experimentalSignal))
    solutionNorms.append(np.linalg.norm(estimation))
plt.scatter(residualNorms, solutionNorms)
plt.xscale("log")
plt.yscale("log")
plt.show()

#Solving the problem with Tikhonov regularization
lam = 0.03 # 0.03
L = np.identity(n)
A = np.concatenate((dict, lam*L))
b = np.concatenate((experimentalSignal, np.zeros(n)))
estimation, r = scipy.optimize.nnls(A, b)
estimatedSignal = np.matmul(dict, estimation)
estimation = estimation[:len(radiuses)] #removing free diffusions
estimation /= radiuses**2 #IMPORTANT STEP
estimation /= np.trapz(estimation, radiuses) #normalization

#quantiles:
def F_estim(r):
    for i in range(len(radiuses)):
        if radiuses[i] > r:
            break
    newY =  estimation[i] + (estimation[i] - estimation[i-1])/(radiuses[i]-radiuses[i-1])
    return scipy.integrate.trapezoid(np.array(list(estimation[:i]) + [newY]), np.array(list(radiuses[:i]) + [float(r)]))

def F_histo(r):
    if r <= histologyRadiuses[0]:
        return 0
    else:
        for i in range(len(histologyRadiuses)):
            if histologyRadiuses[i] > r:
                break
        newY =  histologyPSD[i] + (histologyPSD[i] - histologyPSD[i-1])/(histologyRadiuses[i]-histologyRadiuses[i-1])
        return scipy.integrate.trapezoid(list(histologyPSD[:i]) + [newY], list(histologyRadiuses[:i]) + [float(r)])

histo_q25 = scipy.optimize.fsolve(lambda r : F_histo(r) - 0.25, 10)
histo_q50 = scipy.optimize.fsolve(lambda r : F_histo(r) - 0.5, 10)
histo_q75 = scipy.optimize.fsolve(lambda r : F_histo(r) - 0.75, 20)

estim_q25 = scipy.optimize.fsolve(lambda r : F_estim(r) - 0.25, 10)
estim_q50 = scipy.optimize.fsolve(lambda r : F_estim(r) - 0.5, 10)
estim_q75 = scipy.optimize.fsolve(lambda r : F_estim(r) - 0.75, 20)

print("First quartile:", np.abs(histo_q25 - estim_q25)/histo_q25)
print("Median:", np.abs(histo_q50 - estim_q50)/histo_q50)
print("Third quartile:", np.abs(histo_q75 - estim_q75)/histo_q75)

#plots
fig, axs = plt.subplots(2)
fig.suptitle("Estimation of the Cylinder Radius Distribution by Simulation")

#signals
qAxis = ["0"]
for q in range(len(qValues)):
    qAxis.extend([str(int(1e3 * qValues[q])) + "; 0", str(int(1e3 * qValues[q])) + "; 180"])
axs[0].scatter(qAxis, experimentalSignal, color="b", marker="+")
axs[0].scatter(qAxis, estimatedSignal, color="r", marker="x")
axs[0].legend(["Experimental", "Estimation"])
axs[0].set_xlabel("q=gamma/2pi*G*delta [mm-1]; psi [deg]")
axs[0].set_ylabel("Signal attenuation")
axs[0].grid()

#PSDs
axs[1].plot(histologyRadiuses, histologyPSD, color="b")
axs[1].plot(radiuses, estimation, color="r")
axs[1].legend(["Histology", "Estimation"])
axs[1].set_xlabel("Cylinder radius [um]")
axs[1].set_ylabel("Cylinder radius distribution")
axs[1].grid()

plt.show()

"""
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(results.startingPositions[:, 0], results.startingPositions[:, 1], results.startingPositions[:, 2])
#ax.scatter(results.getEndPositions()[:, 0], results.getEndPositions()[:, 1], results.getEndPositions()[:, 2])
#ax.scatter(results.endPositionsVoxel[:, 0], results.endPositionsVoxel[:, 1], results.endPositionsVoxel[:, 2])
plt.show()
"""
