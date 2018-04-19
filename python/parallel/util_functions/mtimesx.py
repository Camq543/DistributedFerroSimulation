import numpy as np

def mul(arr1,arr2,mode ='N'):
	if mode == 'N':
		return np.einsum('...ij,...jk->...ik',arr1,arr2)
	elif mode == 'T':
		return np.einsum('...ij,...kj->...ik',arr1,arr2)
	elif mode == 'F':
		return np.einsum('...ji,...jk->...ik',arr1,arr2)