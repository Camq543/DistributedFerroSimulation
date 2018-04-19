import numpy as np
import math
from constants import Constants
from util_functions.matToFile import write_csv
import util_functions.voigtTransforms as vt
import util_functions.rotations as rotations
import util_functions.mtimesx as mtimesx

class MonoCrystal:
    __slots__ = ['nbdom', 'AS', 'Pol0', 'dir100', 'dir110', 'dir111', 'Cmonoc3333','dmonoc333','kappa33','epsferro33']

    def __init__(self):
        self.nbdom = Constants.nbdom #double
        self.AS = Constants.AS #double
        self.set_dir()
        R = self.rHelper()
        self.Cmonoc3333 = self.cHelper(R)
        self.dmonoc333 = self.dHelper(R)
        self.kappa33 = self.kHelper(R)
        self.epsferro33 = self.eHelper()
        print('monocrystal setup finished')

    def rHelper(self):
        #finds the rotation matrix R, which is used to calculate the monocrystal tensors
        R = np.zeros((1,self.nbdom,3,3))
        dirR = np.zeros((3,1))
        matCross = np.zeros((3,3))
        toPut = np.zeros((3,3))
        dirRef = np.array([0,0,1]).reshape(3,1)
        for i in range(self.nbdom):
            cosT = np.dot(dirRef.transpose(),self.dir100[0,i])
            dirR = np.array([self.dir100[0,i,1,0],-1*self.dir100[0,i,0,0],0]).reshape(3,1)
            if np.linalg.norm(dirR) != 0:
                dirR = dirR / np.linalg.norm(dirR)

            matCross[0] = np.array([0,0,dirR[1,0]])
            matCross[1] = np.array([0,0,-dirR[0,0]])
            matCross[2] = np.array([-dirR[1,0],dirR[0,0],0])
        
            toPut = np.identity(3) * cosT + (np.outer(dirR,dirR)*(1-cosT)) - (matCross * np.sqrt(1-(cosT*cosT)))
            R[0,i] = toPut
        return R      


    def cHelper(self,R):
        c11 = Constants.C11
        c12 = Constants.C12
        c44 = Constants.C44
        cmonoc = np.array([
            [c11,c12,c12,0,0,0],
            [c12,c11,c12,0,0,0],
            [c12,c12,c11,0,0,0],
            [0,0,0,c44,0,0],
            [0,0,0,0,c44,0],
            [0,0,0,0,0,c44]])
        cmonoc = vt.v66_3333(cmonoc)
        tiled = np.tile(cmonoc,(1,self.nbdom,1,1,1,1))
        return rotations.order_4(tiled,R)

    def dHelper(self,R):
        d31 = Constants.D31
        d33 = Constants.D33
        d15 = Constants.D15
        dmonoc = np.array([
            [0,0,0,0,d15,0],
            [0,0,0,d15,0,0],
            [d31,d31,d33,0,0,0]])
        dmonoc = vt.v36_333(dmonoc)
        tiled = np.tile(dmonoc,(1,self.nbdom,1,1,1))
        return rotations.order_3(tiled,R)

    def kHelper(self,R):
        k33 = Constants.k33
        k11 = Constants.k11
        kappa = np.array([
            [k11,0,0],
            [0,k11,0],
            [0,0,k33]])
        tiled = np.tile(kappa,(1,self.nbdom,1,1))
        return mtimesx.mul(mtimesx.mul(R,tiled,'N'),R,'T')

    def eHelper(self):
        epz0 = Constants.epz0
        matID = np.tile(np.identity(3),(1,self.nbdom,1,1))
        return epz0/2 * (mtimesx.mul(self.dir100,self.dir100,'T') * 3 - matID)

    def set_dir(self):
        self.dir100 = np.zeros((1,6,3,1)) #initialize dir100 array
        self.dir110 = np.zeros((1,12,3,1)) #initialize dir110 array
        self.dir111 = np.zeros((1,8,3,1)) #initialize dir111 array
        #set dir100 values
        self.dir100[0,0] = np.array([[1],[0],[0]])
        self.dir100[0,1] = np.array([[-1],[0],[0]])
        self.dir100[0,2] = np.array([[0],[1],[0]])
        self.dir100[0,3] = np.array([[0],[-1],[0]])
        self.dir100[0,4] = np.array([[0],[0],[1]])
        self.dir100[0,5] = np.array([[0],[0],[-1]])
        #set dir110 values
        self.dir110[0,0] = np.array([[1],[1],[0]])
        self.dir110[0,1] = np.array([[1],[-1],[0]])
        self.dir110[0,2] = np.array([[-1],[1],[0]])
        self.dir110[0,3] = np.array([[-1],[-1],[0]])
        self.dir110[0,4] = np.array([[0],[1],[1]])
        self.dir110[0,5] = np.array([[0],[1],[-1]])
        self.dir110[0,6] = np.array([[0],[-1],[1]])
        self.dir110[0,7] = np.array([[0],[-1],[-1]])
        self.dir110[0,8] = np.array([[1],[0],[1]])
        self.dir110[0,9] = np.array([[1],[0],[-1]])
        self.dir110[0,10] = np.array([[-1],[0],[1]])
        self.dir110[0,11] = np.array([[-1],[0],[-1]])
        self.dir110 *= (1/math.sqrt(2))
        #set dir111 values
        self.dir111[0,0] = np.array([[1],[1],[1]])
        self.dir111[0,1] = np.array([[1],[1],[-1]])
        self.dir111[0,2] = np.array([[1],[-1],[1]])
        self.dir111[0,3] = np.array([[1],[-1],[-1]])
        self.dir111[0,4] = np.array([[-1],[1],[1]])
        self.dir111[0,5] = np.array([[-1],[1],[-1]])
        self.dir111[0,6] = np.array([[-1],[-1],[1]])
        self.dir111[0,7] = np.array([[-1],[-1],[-1]])
        self.dir111 *= (1/math.sqrt(3))

    def data_write(self):
        print('Writing monocrystal data to CSV')
        workDir = "/home/lapis/Documents/FerroProject/python/DataFiles/mono/"
        write_csv(workDir + 'Cmonoc3333',self.Cmonoc3333)
        write_csv(workDir + 'dir100',self.dir100)
        write_csv(workDir + 'dir110',self.dir110)
        write_csv(workDir + 'dir111',self.dir111)
        write_csv(workDir + 'dmonoc333',self.dmonoc333)
        write_csv(workDir + 'epsferro33',self.epsferro33)
        write_csv(workDir + 'kappa33',self.kappa33)
        print('Done')

    