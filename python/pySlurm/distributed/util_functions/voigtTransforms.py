import numpy as np 
import math
#Collection of transformations of symmetric tensors into matrices using Voigt notation
#Also called Mandel notation, Mandel-voigt notation, Nye notation, or kelvin notation

#store matrix to help transforms that is common to all transforms used
transMat = np.array([
	[0,1,2,1,2,0],
	[0,1,2,2,0,1],
	[1,1,1,math.sqrt(2),math.sqrt(2),math.sqrt(2)]])

def v66_3333(C66):
	#go from 6x6 matrix to 3x3x3x3 tensor
	C3333 = np.zeros((3,3,3,3))
	for t in range(6):
		for u in range(6):
			toPut = 1/(transMat[2,t] * transMat[2,u]) * C66[t,u]
			C3333[int(transMat[0,t]),int(transMat[1,t]),int(transMat[0,u]),int(transMat[1,u])] = toPut
			C3333[int(transMat[1,t]),int(transMat[0,t]),int(transMat[0,u]),int(transMat[1,u])] = toPut
			C3333[int(transMat[0,t]),int(transMat[1,t]),int(transMat[1,u]),int(transMat[0,u])] = toPut
			C3333[int(transMat[1,t]),int(transMat[0,t]),int(transMat[1,u]),int(transMat[0,u])] = toPut
	return C3333

def v36_333(d36):
	#go from 3x6 matrix to 3x3x3 tensor
	d333 = np.zeros((3,3,3))
	for i in range(3):
		for j in range(6):
			toPut = 1/(transMat[2,j]) * d36[i,j]
			d333[i,int(transMat[0,j]),int(transMat[1,j])] = toPut
			d333[i,int(transMat[1,j]),int(transMat[0,j])] = toPut
	return d333

def v3333_66(C3333):
	#go from 3x3x3x3 matrix to 6x6 tensor
	#transforms a collection of n x m tensors
	arrSize = C3333.shape
	C66 = np.zeros((arrSize[0],arrSize[1],6,6))
	for n in range(arrSize[0]):
		for m in range(arrSize[1]):
			for i in range(6):
				for j in range(6):
					toPut = transMat[2,i]*transMat[2,j]*C3333[n,m,int(transMat[0,i]),int(transMat[1,i]),int(transMat[0,j]),int(transMat[1,j])]
					C66[n,m,i,j] = toPut
	return C66

def v333_36(d333):
	#go from 3x3x3 tensor to 3x6 matrix
	#transforms a collection of n x m tensors
	arrSize = d333.shape
	d36 = np.zeros((arrSize[0],arrSize[1],3,6))
	for n in range(arrSize[0]):
		for m in range(arrSize[1]):
			for i in range(3):
				for j in range(6):
					toPut = transMat[2,j] * d333[n,m,i,int(transMat[0,j]),int(transMat[1,j])]
					d36[n,m,i,j] = toPut
	return d36

def v33_61(S33):
	#go from 3x3 matrix to 6x1 vector
	#transforms a collection of n x m tensors
	arrSize = S33.shape
	S61 = np.zeros((arrSize[0],arrSize[1],6,1))
	for n in range(arrSize[0]):
		for m in range(arrSize[1]):
			for i in range(6):
				toPut = transMat[2,i] * S33[n,m,int(transMat[0,i]),int(transMat[1,i])]
				S61[n,m,i,0] = toPut
	return S61


	
