import numpy as np
import math

def order_3(d333,R):
	arrSize = d333.shape
	postRot = np.zeros((arrSize[0],arrSize[1],3,3,3))

	for grain in range(arrSize[0]):
		for direc in range(arrSize[1]):
			for i in range(3):
				for j in range(3):
					for k in range(3):
						for m in range(3):
							for n in range(3):
								for o in range(3):
									toPut = postRot[grain,direc,i,j,k] + R[grain,direc,i,m] * R[grain,direc,j,n] * R[grain,direc,k,o] * d333[grain,direc,m,n,o]
									postRot[grain,direc,i,j,k] = toPut
	return postRot

def order_4(C3333,R):
	arrSize = C3333.shape
	postRot = np.zeros((arrSize[0],arrSize[1],3,3,3,3))

	for grain in range(arrSize[0]):
		for direc in range(arrSize[1]):
			for i in range(3):
				for j in range(3):
					for k in range(3):
						for l in range(3):
							for m in range(3):
								for n in range(3):
									for o in range(3):
										for p in range(3):
											toPut = postRot[grain,direc,i,j,k,l] + R[grain,direc,i,m] * R[grain,direc,j,n] * R[grain,direc,k,o] * R[grain,direc,l,p] * C3333[grain,direc,m,n,o,p]
											postRot[grain,direc,i,j,k,l] = toPut
	return postRot

def matrot(angles):
    c1=math.cos(angles[0]);
    s1=math.sin(angles[0]);
    c2=math.cos(angles[1]);
    s2=math.sin(angles[1]);
    c3=math.cos(angles[2]);
    s3=math.sin(angles[2]);

    A_rot = np.zeros((3,3))
    A_rot[0,0] = c1*c3 - s1*c2*s3
    A_rot[0,1] = -c1*s3 - c2*c3*s1
    A_rot[0,2] = s1*s2
    A_rot[1,0] = c3*s1 + c1*c2*s3
    A_rot[1,1] = c1*c2*c3 - s1*s3
    A_rot[1,2] = -c1*s2
    A_rot[2,0] = s2*s3
    A_rot[2,1] = c3*s2
    A_rot[2,2] = c2
    
    return A_rot

