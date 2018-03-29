import numpy as np 
import math

transMat = np.array([
	[0,1,2,1,2,0],
	[0,1,2,2,0,1],
	[1,1,1,math.sqrt(2),math.sqrt(2),math.sqrt(2)]])

def v66_3333(C66):
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
	d333 = np.zeros((3,3,3))
	for i in range(3):
		for j in range(6):
			toPut = 1/(transMat[2,j]) * d36[i,j]
			d333[i,int(transMat[0,j]),int(transMat[1,j])] = toPut
			d333[i,int(transMat[1,j]),int(transMat[0,j])] = toPut
	return d333

def v3333_66(C3333):
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
	arrSize = S33.shape
	S61 = np.zeros((arrSize[0],arrSize[1],6,1))
	for n in range(arrSize[0]):
		for m in range(arrSize[1]):
			for i in range(6):
				toPut = transMat[2,i] * S33[n,m,int(transMat[0,i]),int(transMat[1,i])]
				S61[n,m,i,0] = toPut
	return S61


	
