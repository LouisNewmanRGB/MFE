import numpy as np
import random
import scipy.stats
import scipy.special
import scipy.optimize
import scipy.signal

class Util():

    pickleDir = "numpy saves"

    def getFilePath(func):
        funcName = str(func).split()[1].split(".")[1]
        return "./" + Util.pickleDir + "/" + funcName

    def recursiveMax(listOfLists):
        if type(listOfLists[0]) == list:
            return np.max([Util.recursiveMax(e) for e in listOfLists])
        else:
            return(np.max(np.abs(listOfLists)))

    def getRandomDirection():
        phi = np.random.uniform(0,np.pi*2)
        costheta = np.random.uniform(-1,1)
        theta = np.arccos( costheta )
        return np.array([np.sin( theta) * np.cos( phi ), np.sin( theta) * np.sin( phi ), np.cos( theta )])

    def getRandomU(size):
        return (random.random()-0.5)*size

    def getRandomQuadratic(radius):
        #returns a random number drawn from the interval [0, radius] following a quadradic pdf
        return radius*random.random()**(1/3)

    def getTimeStep(diffusionTime, nStep):
        if nStep == 0:
            print("WARNING: ZERO TIMESTEPS")
            return 1 #irrelevant
        else:
            return diffusionTime / nStep

    def plotPoints(data, plotter, color):
        xArray = np.array([pos[0] for pos in data])
        yArray = np.array([pos[1] for pos in data])
        zArray = np.array([pos[2] for pos in data])
        plotter(xArray, yArray, zArray, color = color)

    def genAlphaList(nTerms, nIter, a):
        """positive roots (solving for alpha) of a*alpha = tan(a*alpha)"""
        alphaList = [(2*k + 1)*np.pi/2 for k in range(nTerms+1)]
        for j in range(nIter):
            for k in range(nTerms+1):
                alphaList[k] = np.arctan(alphaList[k]) + k*np.pi
        return (np.array(alphaList)/a)[1:]

    def genBetaList(nTerms, nIter):
        """positive roots (solving for beta) of beta = tan(beta)"""
        betaList = [(2*k + 1)*np.pi/2 for k in range(nTerms)]
        for j in range(nIter):
            for k in range(nTerms):
                betaList[k] = np.arctan(betaList[k]) + k*np.pi
        betaList[0] = 0
        return np.array(betaList)

    def Jnp(r,n):
        """derivative of the spherical bessel function of order n"""
        return scipy.special.spherical_jn(n, r, derivative=True)

    def Jnp_zeros(n,nt):
        """calculate the zeros of the derivative of the spherical bessel function of order n"""
        zerosj = np.zeros((n+1, nt))
        zerosj[0] = Util.genBetaList(nt+1, nt)[1:]
        points = Util.genBetaList(nt+n, nt)
        points[0] = 0.01
        racines = np.zeros(nt+n)
        for i in range(1,n+1):
            for j in range(nt+n-i):
                foo = scipy.optimize.brentq(Util.Jnp, points[j], points[j+1], (i,))
                racines[j] = foo
            points = racines
            zerosj[i][:nt] = racines[:nt]
        return (zerosj)

    #########################CDFs and PDFs###########################
    def getPDF(D, diffusionTime):
        def f(x):
            return scipy.stats.maxwell.pdf(x, loc = 0, scale = (2*D*diffusionTime)**0.5)
        return f

    def getCDF(D, diffusionTime):
        def F(x):
            return scipy.stats.maxwell.cdf(x, loc = 0, scale = (2*D*diffusionTime)**0.5)
        return F

    def getPDF_sphere(D, diffusionTime, radius, nTerms, nIter):
        alphaList = Util.genAlphaList(nTerms, nIter, radius)

        def f(r):
            if r < radius:
                return 3*r**2/(radius**3) + \
                   2*r*np.sum([np.exp(-D*diffusionTime*alpha**2)*np.sin(alpha*r)*alpha/(np.sin(alpha*radius)**2) for alpha in alphaList])/radius
            return 0
        return f

    def getCDF_sphere(D, diffusionTime, radius, nTerms, nIter):
        alphaList = Util.genAlphaList(nTerms, nIter, radius)

        def F(r):
            return np.heaviside(r - radius, 0) + np.heaviside(radius - r, 1)*( (r/radius)**3 + \
                        2*np.sum([np.exp(-D*diffusionTime*alpha**2)* \
                                  (np.sin(alpha*r) - alpha*r*np.cos(alpha*r))/(alpha*np.sin(alpha*radius)**2) for alpha in alphaList], axis=0)/radius )
        return F

    def getPDF_plane(a):
        def f(x):
            return ( np.heaviside(x+a/2, 1) - np.heaviside(x-a/2, 1) )/a
        return f

    ####################Root mean square displacements#################
    def RMSD_sphere(radius, D, t, nTerms, nIter):
        betaList = Util.Jnp_zeros(1,nTerms)[-1]
        msd = 6*radius**2 / 5 - 12*radius**2 * np.sum([np.exp(-beta**2 *D*t/radius**2) / (beta**4 - 2*beta**2) for beta in betaList])
        return msd**0.5

    ####################Theoretical SGP signals########################
    def getSignal_sphere_inf(radius):
        def s(q):
            qR = q*radius
            return 9*(qR*np.cos(qR) - np.sin(qR))**2 / qR**6
        return s

    def getSignal_sphere_fin(radius, D, t, nNumber, zeroNumber):
        jpZeros = Util.Jnp_zeros(nNumber,zeroNumber)
        def s(q):
            qR = q*radius
            return (3*scipy.special.spherical_jn(1, qR))**2 / qR**2 \
                   + 6*qR**2*np.sum([scipy.special.spherical_jn(n, qR, derivative=True)**2 * \
                    np.sum([(2*n + 1)*alphaNM**2 * np.exp(-alphaNM**2 * D * t / radius**2)/ ( (alphaNM**2 - n**2 - n) * (alphaNM**2 - qR**2)**2)
                            for alphaNM in jpZeros[n] ], axis=0) for n in range(nNumber) ], axis=0)
        return s

    def getSignal_sphere_red(tStar, nNumber, zeroNumber):
        jpZeros = Util.Jnp_zeros(nNumber,zeroNumber)
        def s(qStar):
            qR = np.pi*qStar
            return (3*scipy.special.spherical_jn(1, qR))**2 / qR**2 + 6*qR**2 * np.sum([(scipy.special.spherical_jn(n, qR, derivative=True))**2 * \
                    np.sum([(2*n + 1)*alphaNM**2 * np.exp(-alphaNM**2 * 4*tStar) / ( (alphaNM**2 - n**2 - n)*(alphaNM**2 - qR**2)**2) \
                            for alphaNM in jpZeros[n]], axis=0) for n in range(nNumber)], axis=0)
        return s

    def getSignal_plane_inf(a):
        def s(q):
            qa = q*a
            return 2*(1-np.cos(qa))/qa**2
        return s

    def getSignal_plane_fin(a, D, t, nNumber):
        def s(q):
            qa = q*a
            return 2*(1-np.cos(qa))/qa**2 + 4*qa**2 * \
                   np.sum([(1-(-1)**n * np.cos(qa) ) * np.exp(-(n*np.pi)**2 * D * t/a**2) / ((n*np.pi)**2 - qa**2)**2 for n in range(1,nNumber+1)], axis=0)
        return s

    def getSignal_plane_red(tStar, nNumber):
        def s(qStar):
            qa = 2*np.pi*qStar
            return 2*(1-np.cos(qa))/qa**2 + 4*qa**2 * \
                   np.sum([(1-(-1)**n * np.cos(qa) ) * np.exp(-(n*np.pi)**2 * tStar) / ((n*np.pi)**2 - qa**2)**2 for n in range(1,nNumber+1)], axis=0)
        return s

    def getSignal_cylinder_fin(radius, D, t, nNumber, zeroNumber):
        def s(q):
            qR = q*radius
            kronecker = scipy.signal.unit_impulse(nNumber)
            return (2*scipy.special.jv(1, qR))**2 / qR**2 + 8*qR**2 * np.sum([(scipy.special.jvp(n, qR))**2 / (1 + kronecker[n]) * \
                    np.sum([alphaNM**2 * np.exp(-alphaNM**2 *D*t/radius**2) / ( (alphaNM**2 - n**2) * (alphaNM**2 - qR**2)**2) \
                            for alphaNM in scipy.special.jnp_zeros(n, zeroNumber)], axis=0) for n in range(nNumber)], axis=0)
        return s

    def getSignal_cylinder_red(tStar, nNumber, zeroNumber):
        def s(qStar):
            qR = np.pi*qStar
            kronecker = scipy.signal.unit_impulse(nNumber)
            return (2*scipy.special.jv(1, qR))**2 / qR**2 + 8*qR**2 * np.sum([(scipy.special.jvp(n, qR))**2 / (1 + kronecker[n]) * \
                    np.sum([alphaNM**2 * np.exp(-alphaNM**2 *4*tStar) / ( (alphaNM**2 - n**2) * (alphaNM**2 - qR**2)**2) \
                            for alphaNM in scipy.special.jnp_zeros(n, zeroNumber)], axis=0) for n in range(nNumber)], axis=0)
        return s
