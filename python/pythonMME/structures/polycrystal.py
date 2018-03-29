import numpy as np
from constants import Constants
from util_functions.matToFile import write_csv
import util_functions.rotations as rotations
import util_functions.voigtTransforms as vt
import util_functions.eshelby as es
import util_functions.mtimesx as mtimesx



class PolyCrystal:
    __slots__ = ['nbg','monoc','angles','mRot','AS',
    'dir100','dir110','dir111','Pol','epsferro61','kappa33',
    'd36','C66','invC66','ceffMeca','ceffElec','eshMeca',
    'eshElec','cStarMeca','cStarElec']

    def __init__(self, monoc):
        #import monocrystal data and construct grain orientations
        self.monoc = monoc
        self.angles = self.anglesetter()
        self.nbg = len(self.angles)
        self.AS = monoc.AS

        #using the grain orientations, rotate monocrystal 
        #arrays to be in line with each grain
        self.mRot = np.zeros((self.nbg,1,3,3))
        for i in range(self.nbg):
            self.mRot[i,0] = rotations.matrot(self.angles[i,0])
        self.dir100 = mtimesx.mul(self.mRot,monoc.dir100,'N')
        self.dir110 = mtimesx.mul(self.mRot,monoc.dir110,'N')
        self.dir111 = mtimesx.mul(self.mRot,monoc.dir111,'N')
        self.Pol = self.dir100 * Constants.Pol0
        self.epsferro61 = vt.v33_61(mtimesx.mul(mtimesx.mul(self.mRot,monoc.epsferro33,'N'),self.mRot,'T'))
        self.kappa33 = mtimesx.mul(mtimesx.mul(self.mRot,monoc.kappa33,'N'),self.mRot,'T')
        self.d36 = vt.v333_36(rotations.order_3(np.tile(monoc.dmonoc333,(self.nbg,1,1,1,1)),np.tile(self.mRot,(1,monoc.nbdom,1,1))))
        self.C66 = vt.v3333_66(rotations.order_4(np.tile(monoc.Cmonoc3333,(self.nbg,1,1,1,1,1)),np.tile(self.mRot,(1,monoc.nbdom,1,1))))
        self.invC66 = np.zeros((self.nbg,monoc.nbdom,6,6))
        for n in range(self.nbg):
            for j in range(monoc.nbdom):
                self.invC66[n,j] = np.linalg.inv(self.C66[n,j])

        #setup homoginization tensors
        self.ceffMeca = self.C66[0,0]
        self.ceffElec = self.kappa33[0,0]
        self.eshMeca = es.Eshelby66_muphylin(self.ceffMeca,1)
        self.eshElec = es.Eshelby33_muphylin(self.ceffElec,1)
        self.cStarMeca = np.tile(np.dot(self.ceffMeca,np.linalg.inv(self.eshMeca)) - self.ceffMeca,(self.nbg,1,1,1))
        self.cStarElec = np.tile(np.dot(self.ceffElec,np.linalg.inv(self.eshElec)) - self.ceffElec,(self.nbg,1,1,1))
        self.eshMeca = np.tile(self.eshMeca,(self.nbg,1,1,1))
        self.eshElec = np.tile(self.eshElec,(self.nbg,1,1,1))

        print('Polycrystal setup finished')


    def anglesetter(self):
        n1 = Constants.n1
        n2 = Constants.n2
        n3 = Constants.n3
        toReturn = np.zeros((n1*n2*n3,1,3,1))

        angle1 = np.linspace(0,2*np.pi,n1 + 1)[:-1]
        angle2 = np.arccos(np.linspace(0,1,n2 + 1)[:-1])
        angle3 = np.linspace(0,np.pi/2,n3 + 1)[:-1]
        
        count = 0
        for i in range(n1):
            for j in range(n2):
                for k in range(n3):
                    toReturn[count,0] = np.array([angle1[i],angle2[j],angle3[k]]).reshape(3,1)
                    count += 1
        return toReturn

    def data_write(self):
        print('Writing polycrystal data to CSV')
        workDir = "/home/lapis/Documents/FerroProject/python/DataFiles/poly/"
        write_csv(workDir + 'angles',self.angles)
        write_csv(workDir + 'C66',self.C66)
        write_csv(workDir + 'cStarElec',self.cStarElec)
        write_csv(workDir + 'cStarMeca',self.cStarMeca)
        write_csv(workDir + 'd36',self.d36)
        write_csv(workDir + 'epsferro61',self.epsferro61)
        write_csv(workDir + 'eshElec',self.eshElec)
        write_csv(workDir + 'eshMeca',self.eshMeca)
        write_csv(workDir + 'invC66',self.invC66)
        write_csv(workDir + 'kappa33',self.kappa33)
        write_csv(workDir + 'M_rot',self.mRot)
        print('Done')




