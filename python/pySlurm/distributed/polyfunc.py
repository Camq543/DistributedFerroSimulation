import numpy as np
import util_functions.mtimesx as mtimesx

def polyfunc(polyc, Tmacro, Emacro, Tlocinit=None, Elocinit=None):
    nbg = polyc.nbg

    Smacro = np.zeros((6,1)) #initial response is 0 vector
    Dmacro = np.zeros((3,1))
    if(Tlocinit is None or Elocinit is None):
        #If no tensor is passed in, construct them from starting vectors
        T61 = np.tile(Tmacro,(nbg,1,1,1))
        E = np.tile(Emacro,(nbg,1,1,1))
    else:
        #Else set tensors we use for calculation equal to initial args
        T61 = Tlocinit 
        E = Elocinit

    if (np.linalg.norm(Tmacro) + np.linalg.norm(Emacro) == 0):
        #If initial stress and electric field vectors are 0, no calculation is necessary
        Tloc = np.tile(Tmacro,(nbg,1,1,1))
        Eloc = np.tile(Emacro,(nbg,1,1,1))
        return (Smacro, Dmacro, Tloc, Eloc,0)
    #initialize mechanical and electric convergence variables to infinity
    mecConvergence = float('inf')
    elecConvergence = float('inf')
    convergence = float('inf')
    mecConvMem = float('inf')
    elecConvMem = float('inf')
    
    #alpha list stores a list of increasing values where 0 < x < 1
    #these values determine the degree to which each iteration changes the calculated tensors
    #a higher alpha value means the previous tensors will change less
    #if an iteration fails to move closer to convergence, we increase the value of alpha
    alphaNum = 20
    alphaList = np.logspace(np.log10(.9),-2,alphaNum)
    alphaList = 1 - alphaList
    #initialize indexes for alphaList
    mecIndex = 0
    elecIndex = 0
    alphaMec = alphaList[mecIndex]
    alphaElec = alphaList[elecIndex]
    
    count = 0 #store number of iterations
    while convergence > 1e-3:
        #calculation converges when the change between iterations is less than .1%
        #convergence criteria can be adjusted to change precision
        count +=1 

        #store mechanical and electric response tensor
        S61,D = LDC_monoc_vectorize(polyc,T61,E)
        #hold the mean value for S61 and D to create response vectors
        Smacro = np.mean(S61,0) 
        Dmacro = np.mean(D,0) 

        #hold onto previous tensors in order to calculate convergence
        Tmem = T61
        Emem = E

        #calculate new mechanical and induction tensors
        T61 = np.tile(Tmacro,(nbg,1,1,1)) + (mtimesx.mul(polyc.cStarMeca,np.tile(Smacro,(nbg,1,1,1)) - S61[:,np.newaxis,:,:],'N'))
        E = np.tile(Emacro,(nbg,1,1,1)) + (mtimesx.mul(polyc.cStarElec,np.tile(Dmacro,(nbg,1,1,1)) - D[:,np.newaxis,:,:],'N'))
        
        #each tensor changes by a percent equal to one minus the current alpha value
        T61 = Tmem*alphaMec + T61 * (1-alphaMec)
        E = Emem*alphaElec + E * (1-alphaElec)
        #calculate difference with original tensors
        diffMec = T61 - Tmem
        diffElec = E - Emem
        #hold onto previous convergence values to make sure simulation isnt diverging
        mecConvMem = mecConvergence
        elecConvMem = elecConvergence
        #set current convergence based on percent change in the norm of the tensors
        mecConvergence = np.linalg.norm(diffMec)/(1 + np.linalg.norm(T61))
        elecConvergence = np.linalg.norm(diffElec)/(1 + np.linalg.norm(E))

        if mecConvergence > mecConvMem:
            #if mechanic response is not converging, increase alpha for mechanics
            mecIndex += 1
            if mecIndex < alphaNum:
                alphaMec = alphaList[mecIndex]
            else:
                #if we reached max alpha value, simulation failed to converge
                print("No convergence for mechanics")
                break

        if elecConvergence > elecConvMem:
            #if electric response is not converging, increase alpha for induction
            elecIndex += 1
            if elecIndex < alphaNum:
                alphaElec = alphaList[elecIndex]
            else:
                #if we reached max alpha value, simulation failed to converge
                print("No convergence for electricity")
                break

        convergence = max(mecConvergence,elecConvergence) #set convergence check value to whichever tensor converged less
        # if count >=400:
        #     #if iteration count is too high, end simulation. Criteria can be changed to adjust precision and computation speed
        #     print("Convergence ERROR, max count reached")
        #     print("Tmacro", Tmacro)
        #     print("Emacro", Emacro)
        #     convergence = 0

    #set return values
    Tloc = T61
    Eloc = E
    Smacro = np.mean(S61,0)
    Dmacro = np.mean(D,0)

    return (Smacro, Dmacro, Tloc, Eloc, count) #return respose vectors and final calculated tensors    


       

def LDC_monoc_vectorize(polyc, T61, E):
    #based on input poly crystal and stress/electric tensor, determine response vectors
    #I don't really know how to explain the physics, but this function is basically a bunch of matrix multiplication
    #this is where the vast majority of computation time is spent
    dpz = mtimesx.mul(polyc.d36,T61,'N')
    epspz = mtimesx.mul(polyc.d36,E,'F')
    
    W = -mtimesx.mul(E,polyc.Pol,'F') - mtimesx.mul(T61,polyc.epsferro61,'F') - 2*mtimesx.mul(E,dpz,'F')    
    
    fvol = np.exp(-polyc.AS * W) / np.tile(np.sum(np.exp(-polyc.AS * W),1)[:,np.newaxis,:,:],(1,6,1,1))
    
    D = np.sum(mtimesx.mul((mtimesx.mul(polyc.kappa33,E,'N') + dpz + polyc.Pol),fvol,'N'),1)

    S61 = np.sum(mtimesx.mul((mtimesx.mul(polyc.invC66,T61,'N') + epspz + polyc.epsferro61),fvol,'N'),1)

    Wtot = np.sum(np.sum(mtimesx.mul(W,fvol,'N'),1),0)

    return(S61,D)