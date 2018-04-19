import numpy as np
#Simple implementation of mtimesx library for matlab (now called ND array multiplication I think)

def mul(arr1,arr2,mode ='N'):
	#Given two input tensors of size (M x N x A) and (I x J x B), where A and B are matrices
	#Calculate output tensor of size (max(M,I) x max(N,J) x C), 
	#where C is the product of matrix multiplication of A and B for each matrix in the tensor
	#Does not compute all possible combinations of matrix multiplications as in tensor multiplication,
	#but instead computes product for each matrix applying to the same value of M and I and of N and J
	#In the event M is not equal to I, or N is not equal to J, broadcasts the smaller of the two 
	#to create a final tensor of size (max(M,I) x max(N,J) x C)
	if mode == 'N':
		#mode of N leaves each individual matrix as is and computes product so C = A x B
		return np.einsum('...ij,...jk->...ik',arr1,arr2)
	elif mode == 'T':
		#mode of T transposes each matrix in the second tensor so C = A x B^t
		return np.einsum('...ij,...kj->...ik',arr1,arr2)
	elif mode == 'F':
		#mode of F transposes each matrix in the first tensor so C = A^t x B
		return np.einsum('...ji,...jk->...ik',arr1,arr2)