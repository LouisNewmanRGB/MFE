import numpy as np
import matplotlib.pyplot as plt

def plotMaskedImages(mask, images, plotTitle):
    if plotTitle != None:
        masked = np.ma.masked_where(mask == 0, mask)
        for b in bValues[1:]:
            for psi in [0, 180]:
                plt.imshow(np.transpose(images[(b, psi)][:, :, 2]), interpolation='none', cmap='gray')
                plt.imshow(np.transpose(masked), cmap='jet', interpolation='none', alpha=0.5)
                plt.title(plotTitle)
                #plt.title(plotTitle + "\nb={b}, psi={psi}".format(b=b, psi=psi))
                plt.show()

def createImageDict(dir, bValues, plotImages=False):
    images = {}
    images[0] = np.load(dir + "0_1_mean.npy")
    for b in bValues[1:]:
        images[(b, 0)] = (np.load(dir + "{b}_1_mean.npy".format(b=b)) + np.load(dir + "{b}_3_mean.npy".format(b=b)))/2
        images[(b, 180)] = (np.load(dir + "{b}_2_mean.npy".format(b=b)) + np.load(dir + "{b}_4_mean.npy".format(b=b)))/2
    return images

def createIntensityMask(slice, images, intensityThreshold):
    intensity = images[0].copy()
    intensity = intensity[:, :, slice]
    mask = intensity > intensityThreshold*np.max(intensity)
    return mask

def createMask(slice, images, diffbValue, intensityThreshold, diffThreshold):
    img0 = images[(diffbValue, 0)].copy()
    img180 = images[(diffbValue, 180)].copy()
    img0 = img0[:, :, slice]
    img180 = img180[:, :, slice]
    diff = np.abs(img0 - img180)

    intensity = createIntensityMask(slice, images, intensityThreshold)
    diff = diff > diffThreshold*np.max(diff)
    mask = intensity*diff
    return mask

def createMaskFA(slice, dir, images, intensityThreshold, FAthreshold):
    FA = np.load(dir + "fa_affine_registered.npy")[:, :, slice]/1000
    FA_x = np.load(dir + "color_fa_x_affine_registered.npy")[:, :, slice]
    FA_y = np.load(dir + "color_fa_y_affine_registered.npy")[:, :, slice]
    FA_z = np.load(dir + "color_fa_z_affine_registered.npy")[:, :, slice]

    mask = createIntensityMask(slice, images, intensityThreshold)
    mask *= FA_x < FA_z
    mask *= FA_y < FA_z
    mask *= FA > FAthreshold

    return mask

def getSignal(slice, images, mask, bValues):
    signalValues = np.empty(2*(len(bValues) - 1) + 1)
    signalValues[0] = np.mean(images[0][:, :, slice][np.where(mask != 0)])
    for b in range(1, len(bValues)):
        bValue = bValues[b]
        m1 = np.mean(images[bValue, 0][:, :, slice][np.where(mask != 0)])
        m2 = np.mean(images[bValue, 180][:, :, slice][np.where(mask != 0)])
        signalValues[2*b - 1] = m2
        signalValues[2*b] = m1
    signalValues *= 1000/signalValues[0]
    return signalValues

def getAverageSignalDiff(dir, bValues, intensityThreshold=0.4, diffThreshold=0.4, plotTitle=None):
    images = createImageDict(dir, bValues)
    signal = np.empty((5, 2*len(bValues)-1))
    for i in range(5):
        mask = createMask(i, images, 3500, intensityThreshold, diffThreshold)
        signal[i,:] = getSignal(i, images, mask, bValues)
        if i == 2:
            plotMaskedImages(mask, images, plotTitle)
    return np.mean(signal, axis=0), np.std(signal, axis=0)

def getAverageSignalFA(dir, bValues, intensityThreshold=0.3, FAthreshold=0.275, plotTitle=None):
    images = createImageDict(dir, bValues)
    signal = np.empty((5, 2*len(bValues)-1))
    for i in range(5):
        mask = createMaskFA(i, dir, images, intensityThreshold, FAthreshold)
        signal[i,:] = getSignal(i, images, mask, bValues)
        if i == 2:
            plotMaskedImages(mask, images, plotTitle)
    return np.mean(signal, axis=0), np.std(signal, axis=0)

def getXaxisValues(bValues, delta=12, Delta=60):
    qValues = np.sqrt(1e3*bValues/(8*np.pi**2 *(Delta - delta/3)))
    xAxisValues = [str(0)]
    for b in range(1, len(bValues)):
        qValue = qValues[b]
        xAxisValues.extend([str(int(np.round(qValue))) + "; 0", str(int(np.round(qValue))) + "; 180"])
    return xAxisValues

#17 march
bValues = np.array([0, 220, 870, 2000, 3500, 5400])
xAxisValues = getXaxisValues(bValues)

#1
dir1 = "mri images/delhaize/"
signals1diff, stds1diff = getAverageSignalDiff(dir1, bValues)#, plotTitle="1 Delhaize Diff")
signals1FA, stds1FA = getAverageSignalFA(dir1, bValues, plotTitle="MRI Scan of Asparagus Using the DDE Sequence"
                                                                  "\nFractional Anisotropy-based Mask"
                                                                  "\nb=0.22ms/um2 psi=0deg")
                                         #plotTitle="1 Delhaize FA")
#signals1FA, stds1FA = getAverageSignalFA(dir1, bValues)

#2
dir2 = "mri images/okay/"
signals2diff, stds2diff = getAverageSignalDiff(dir2, bValues)
signals2FA, stds2FA = getAverageSignalFA(dir2, bValues, plotTitle="2 Okay")
#signals2FA, stds2FA = getAverageSignalFA(dir2, bValues)

#export
np.save("numpy saves/asparagus_signals.npy", np.array([signals1diff, signals1FA, signals2diff, signals2FA]))

#plot
ducheneSignal = 1000* np.array([1, 0.8, 0.84, 0.4, 0.51, 0.144, 0.248, 0.048, 0.112, 0.02, 0.056])
#plt.errorbar(xAxisValues, signals1diff, stds1diff, fmt="_")
#plt.errorbar(xAxisValues, signals1FA, stds1FA, fmt="_")
#plt.errorbar(xAxisValues, signals2diff, stds2diff, fmt="_")
#plt.errorbar(xAxisValues, signals2FA, stds2FA, fmt="_")
#plt.scatter(xAxisValues, ducheneSignal)
plt.scatter(xAxisValues, 1 - signals1diff/ducheneSignal, marker="+")
plt.scatter(xAxisValues, 1 - signals1FA/ducheneSignal, marker="x")
plt.scatter(xAxisValues, 1 - signals2diff/ducheneSignal, marker="+")
plt.scatter(xAxisValues, 1 - signals2FA/ducheneSignal, marker="x")
plt.legend(["Pack 1 - Difference", "Pack 1 - FA", "Pack 2 Difference", "Pack 2 - FA"])
plt.xlabel("q=gamma*G*delta [mm-1]; psi [deg]")
#plt.ylabel("MRI signal attenuation")
plt.ylabel("(Sduc-S)/Sduc")
#plt.title("Experimental application of DDE sequence on asparaguses")
plt.title("Relative difference between our experimental MRI signals and Duchene's")
plt.grid()
plt.show()

#DDE gradient directions
#+1 +1 0 +1 +1 0
#+1 +1 0 -1 -1 0
#+1 -1 0 +1 -1 0
#+1 -1 0 -1 +1 0)
