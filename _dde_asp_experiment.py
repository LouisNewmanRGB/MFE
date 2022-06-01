import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.optimize
import scipy.stats

from SimulationCpp import Simulation
from DDEasp import DDEasp
from SimulationCppResults import SimulationCppResults

simulateDict = False
nStep = 50
nPart = int(5*1e5)
delta = 12 #ms
Delta = 60 #60 #ms
diffusionTime = 2*Delta + 2*delta
T2 = np.inf
D = 2 #um2/ms

qValues = 1e-3 * np.array([7, 14, 21, 28, 35])#[:-2] #mm-1
angles = [0, 180]

#data
data = np.load("numpy saves/asparagus_signals.npy")/1000
signals1diff = data[0]#[:-4]
signals1FA = data[1]#[:-4]
signals2diff = data[2]#[:-4]
signals2FA = data[3]#[:-4]
experimentalSignals = [signals1diff, signals1FA, signals2diff, signals2FA]

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
    #np.save("numpy saves/AMICO_dict_experiment.npy", dict)
else:
    dict = np.load("numpy saves/AMICO_dict_duchene_500k.npy")#[:-4, :]

#optimization: solving NNLS problem
#lambdas found:
chosenLambdas = [0.03, 0.03, 0.03, 0.03]

estimations = np.empty((4, len(radiuses)))
signalEstimations = np.empty((4, 2*len(qValues)+1))
for k in range(4):
    experimentalSignal = experimentalSignals[k]
    n = len(radiuses) + len(diffusionCoefs)

    #L-curve method -> find the "corner"
    lambdas = [0.001, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05, 0.06, 0.1]#, 0.2, 0.5, 0.7, 0.9, 1, 2]
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
    plt.xlabel("Residual Norm")
    plt.ylabel("Solution Norm")
    plt.xscale("log")
    plt.yscale("log")
    plt.show()

    #Solving the problem with Tikhonov regularization
    lam = chosenLambdas[k]
    L = np.identity(n)
    A = np.concatenate((dict, lam*L))
    b = np.concatenate((experimentalSignal, np.zeros(n)))
    estimation, r = scipy.optimize.nnls(A, b)
    estimatedSignal = np.matmul(dict, estimation)
    estimation = estimation[:len(radiuses)] #removing free diffusions
    estimation /= radiuses**2 #IMPORTANT STEP
    estimation /= np.trapz(estimation, radiuses) #normalization

    estimations[k] = estimation
    signalEstimations[k] = estimatedSignal

#plots
fig, axs = plt.subplots(2)
fig.suptitle("Estimation of the Cylinder Size Distribution by Simulation")
markers = ["+", "x"]
for k in range(4):
    #signals
    qAxis = ["0"]
    for q in range(len(qValues)):
        qAxis.extend([str(int(1e3 * qValues[q])) + "; 0", str(int(1e3 * qValues[q])) + "; 180"])
    m = markers[k % 2]
    axs[0].scatter(qAxis, experimentalSignals[k], marker=m)

    #PSDs
    axs[1].plot(radiuses, estimations[k])

axs[0].legend(["Pack 1 - Difference", "Pack 1 - FA", "Pack 2 Difference", "Pack 2 - FA"])
axs[0].set_xlabel("q=gamma*G*delta [mm-1]; psi [deg]")
axs[0].set_ylabel("Signal attenuation")
axs[0].grid()

#Histology data
names = ["D3B (2)", "O2H (2)"]
#names = ["O2H (2)"]
for name in names:
    df = pd.read_excel("histology/{name}.xlsx".format(name=name))
    df = np.array(df)
    histologyRadiuses = df[:,-1]
    #histologyRadiuses *= 1.11 #volume reduction correction?
    bw = 2 *      2*scipy.stats.iqr(histologyRadiuses, rng=(25, 75))/(len(histologyRadiuses))**(1/3)
    nBins = int((np.max(histologyRadiuses) - np.min(histologyRadiuses))/bw)
    histY, histX = np.histogram(histologyRadiuses, bins=nBins, density = True)
    histX = (histX[:-1] + histX[1:])/2
    axs[1].plot(histX, histY, linestyle='dashed')

#quantiles:
def F_estim(r, estimation):
    for i in range(len(radiuses)):
        if radiuses[i] > r:
            break
    newY =  estimation[i] + (estimation[i] - estimation[i-1])/(radiuses[i]-radiuses[i-1])
    return scipy.integrate.trapezoid(np.array(list(estimation[:i]) + [newY]), np.array(list(radiuses[:i]) + [float(r)]))

def F_histo(r):
    if r <= histX[0]:
        return 0
    else:
        for i in range(len(histX)):
            if histX[i] > r:
                break
        newY =  histY[i] + (histY[i] - histY[i-1])/(histX[i]-histX[i-1])
        return scipy.integrate.trapezoid(list(histY[:i]) + [newY], list(histX[:i]) + [float(r)])

error25, error50, error75 = [], [], []
for k in range(4):
    histo_q25 = scipy.optimize.fsolve(lambda r : F_histo(r) - 0.25, 10)
    histo_q50 = scipy.optimize.fsolve(lambda r : F_histo(r) - 0.5, 10)
    histo_q75 = scipy.optimize.fsolve(lambda r : F_histo(r) - 0.75, 20)

    estim_q25 = scipy.optimize.fsolve(lambda r : F_estim(r, estimations[k]) - 0.25, 10)
    estim_q50 = scipy.optimize.fsolve(lambda r : F_estim(r, estimations[k]) - 0.5, 10)
    estim_q75 = scipy.optimize.fsolve(lambda r : F_estim(r, estimations[k]) - 0.75, 20)

    error25.append(np.abs(histo_q25 - estim_q25)/histo_q25)
    error50.append(np.abs(histo_q50 - estim_q50)/histo_q50)
    error75.append(np.abs(histo_q75 - estim_q75)/histo_q75)

print("First quartile:", error25, np.mean(error25))
print("Median:", error50, np.mean(error50))
print("Third quartile:", error75, np.mean(error75))

axs[1].legend(["Pack 1 - Difference", "Pack 1 - FA", "Pack 2 Difference", "Pack 2 - FA"] + ["Histology - asparagus 1", "Histology - asparagus 2"]) #["Histology " + name for name in names])
axs[1].set_xlabel("Cylinder radius [um]")
axs[1].set_ylabel("Cylinder radius distribution")
axs[1].grid()

plt.show()
