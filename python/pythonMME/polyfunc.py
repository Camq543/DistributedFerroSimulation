import numpy as np
import util_functions.mtimesx as mtimesx

def polyfunc(polyc, Tmacro, Emacro, Tlocinit=None, Elocinit=None):

    nbg = polyc.nbg

    Smacro = np.zeros((6,1))
    Dmacro = np.zeros((3,1))
    Tloc = np.zeros((nbg,1,6,1))
    Eloc = np.zeros((nbg,1,3,1))
    T61 = Tlocinit
    E = Elocinit

    if (np.linalg.norm(Tmacro) + np.linalg.norm(Emacro) == 0):
        Tloc = np.tile(Tmacro,(nbg,1,1,1))
        Eloc = np.tile(Emacro,(nbg,1,1,1))
        return (Smacro, Dmacro, Tloc, Eloc,0)
    else:
        Tloc = Tlocinit
        Eloc = Elocinit

    mecConvergence = float('inf')
    elecConvergence = float('inf')
    convergence = float('inf')
    mecConvMem = float('inf')
    elecConvMem = float('inf')
    alphaNum = 20
    alphaList = np.logspace(np.log10(.9),-2,alphaNum)
    alphaList = 1 - alphaList

    mecIndex = 0
    elecIndex = 0
    alphaMec = alphaList[mecIndex]
    alphaElec = alphaList[elecIndex]
    count = 0

    while convergence > 1e-3:
        count +=1

        S61,D = LDC_monoc_vectorize(polyc,T61,E)
        Smacro = np.mean(S61,0)
        Dmacro = np.mean(D,0)

        Tmem = T61
        Emem = E
        T61 = np.tile(Tmacro,(nbg,1,1,1)) + (mtimesx.mul(polyc.cStarMeca,np.tile(Smacro,(nbg,1,1,1)) - S61[:,np.newaxis,:,:],'N'))
        E = np.tile(Emacro,(nbg,1,1,1)) + (mtimesx.mul(polyc.cStarElec,np.tile(Dmacro,(nbg,1,1,1)) - D[:,np.newaxis,:,:],'N'))
        
        T61 = Tmem*alphaMec + T61 * (1-alphaMec)
        E = Emem*alphaElec + E * (1-alphaElec)

        diffMec = T61 - Tmem
        diffElec = E - Emem

        mecConvMem = mecConvergence
        elecConvMem = elecConvergence

        mecConvergence = np.linalg.norm(diffMec)/(1 + np.linalg.norm(T61))
        elecConvergence = np.linalg.norm(diffElec)/(1 + np.linalg.norm(E))

        if mecConvergence > mecConvMem:
            mecIndex += 1
            if mecIndex < alphaNum:
                alphaMec = alphaList[mecIndex]
            else:
                print("No convergence for mechanics")
                count = 400

        if elecConvergence > elecConvMem:
            elecIndex += 1
            if elecIndex < alphaNum:
                alphaElec = alphaList[elecIndex]
            else:
                print("No convergence for electricity")
                count = 400

        convergence = max(mecConvergence,elecConvergence)
        #if count >=400:
        #    print("Convergence ERROR, max count reached")
        #    print("Tmacro", Tmacro)
        #    print("Emacro", Emacro)
        #    convergence = 0

    Tloc = T61
    Eloc = E
    Smacro = np.mean(S61,0)
    Dmacro = np.mean(D,0)

    #print("Convergence")
    return (Smacro, Dmacro, Tloc, Eloc, count)    


       

def LDC_monoc_vectorize(polyc, T61, E):
    dpz = mtimesx.mul(polyc.d36,T61,'N')
    epspz = mtimesx.mul(polyc.d36,E,'F')
    
    W = -mtimesx.mul(E,polyc.Pol,'F') - mtimesx.mul(T61,polyc.epsferro61,'F') - 2*mtimesx.mul(E,dpz,'F')    
    
    fvol = np.exp(-polyc.AS * W) / np.tile(np.sum(np.exp(-polyc.AS * W),1)[:,np.newaxis,:,:],(1,6,1,1))
    
    #D = np.tile(np.sum(mtimesx.mul((mtimesx.mul(polyc.kappa33,E,'N') + dpz + polyc.Pol),fvol,'N'),1)[:,np.newaxis,:,:],(1,6,1,1))
    D = np.sum(mtimesx.mul((mtimesx.mul(polyc.kappa33,E,'N') + dpz + polyc.Pol),fvol,'N'),1)

    #S61 = np.tile(np.sum(mtimesx.mul((mtimesx.mul(polyc.invC66,T61,'N') + epspz + polyc.epsferro61),fvol,'N'),1)[:,np.newaxis,:,:],(1,6,1,1))
    S61 = np.sum(mtimesx.mul((mtimesx.mul(polyc.invC66,T61,'N') + epspz + polyc.epsferro61),fvol,'N'),1)

    Wtot = np.sum(np.sum(mtimesx.mul(W,fvol,'N'),1),0)

    return(S61,D)