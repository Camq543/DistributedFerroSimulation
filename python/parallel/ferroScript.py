#python/3rd party libs
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing
import time
from functools import partial

#my functions
from util_functions.matToFile import write_csv
from constants import Constants
from polyfunc import polyfunc
from structures.monocrystal import MonoCrystal as Monoc
from structures.polycrystal import PolyCrystal as Polyc

def main():
	workDir = "/home/lapis/Documents/FerroProject/python/pythonMME/DataFiles/"
	mono = Monoc()
#	mono.data_write()
	poly = Polyc(mono)
#	poly.data_write()

	n = 51
	#n = 1
	Evalue = np.linspace(0,2e6,n)
	Tvalue = np.array([-100,-50,0,50,100]) * 1e6
	#Evalue = np.array([1e6])
	#Tvalue = np.array([10e7])
	m = len(Tvalue)

	t1 = time.time()
	pool = multiprocessing.Pool()
	func = partial(run_poly,poly,Tvalue,Evalue)
	results = pool.map(func, range(m))

	polycD = np.zeros((m,n))
	polycS = np.zeros((m,n))

	for k in range(m):
		polycD[k,:] = results[k][1][:]
		polycS[k,:] = results[k][2][:]

	pool.close()
	polyTime = time.time() - t1
	print("Polycrystal function took " + str(polyTime) + " seconds.")
	# outFile.close()

	fig, ax = plt.subplots()
	colors = ['c','r','b','g','m']
	for k in range(m):
		xvals = []
		yvals = []
		for i in range(n):
			xvals.append(Evalue[i])
			yvals.append(polycD[k,i])
		ax.plot(xvals,yvals,colors[k])

	fig2,ax2 = plt.subplots()
	for k in range(m):
		xvals = []
		yvals = []
		for i in range(n):
			xvals.append(Evalue[i])
			yvals.append(polycS[k,i] - polycS[k,0])
		ax2.plot(xvals,yvals,colors[k])

	#plt.show()

	nbt = Constants.nbt
	nbe = Constants.nbe
	LApp = np.zeros((nbt,nbe,9,9))
	Tvalue = np.linspace(0,Constants.tSpace,nbt) * 1e6
	Evalue = np.linspace(0,Constants.eSpace,nbe) * 1e6
	data = []

	t1 = time.time()
	for i in range(nbt):
		data.append(map_piezo_e(poly,Tvalue,Evalue,nbe,i))

	# print(len(data))
	# print(len(data[0]))
	# print(len(data[0][0]))
	# print(data[0][0][0][1].shape)
	for tDat in data:
		for eDat in tDat[1]:
			for jDat in eDat[1]:
				LApp[tDat[0],eDat[0],:,jDat[0]] = jDat[1][:,0]


	write_csv("LApp",LApp)

	piezoTime = time.time() - t1
	print("Piezo coefficents took " + str(piezoTime) + " seconds.")

def run_poly(poly,Tvalue,Evalue,k):
	n = len(Evalue)
	tempD = np.zeros(n)
	tempS = np.zeros(n)

	Tmacro = np.array([0,0,Tvalue[k],0,0,0]).reshape(6,1)

	Tloc = np.tile(Tmacro,(poly.nbg,1,1,1))
	Eloc = np.zeros((poly.nbg,1,3,1))

	for i in range(n):
		Emacro = np.array([0,0,Evalue[i]]).reshape(3,1)

		Smacro,Dmacro,Tloc,Eloc,iterCount = polyfunc(poly,Tmacro,Emacro,Tloc,Eloc)
		# outFile.write(str(i) + "," + str(k) + "," + str(iterCount) + "\n")

		tempD[i] = Dmacro[2,:]
		tempS[i] = Smacro[2,:]

	return k,tempD,tempS

def map_piezo_e(poly,Tvalue,Evalue,nbe,i):
	pool = multiprocessing.Pool()
	func = partial(run_piezo, poly, Tvalue, Evalue,i)
	data = pool.map(func,range(nbe))
	pool.close()
	return i,data



def run_piezo(poly,Tvalue,Evalue,i,k):
	alph = Constants.alpha
	mecVar = np.zeros((6,1))
	elecVar = np.zeros((3,1))
	delta = 0
	TStat = np.array([0,0,Tvalue[k],0,0,0]).reshape(6,1)
	EStat = np.array([0,0,Evalue[i]]).reshape(3,1)
	Tloc = np.tile(TStat,(poly.nbg,1,1,1))
	Eloc = np.tile(EStat,(poly.nbg,1,1,1))
	toReturn = []
	# print(TStat)
	# print(EStat)

	Smacro,Dmacro,Tinit,Einit,iterCount = polyfunc(poly,TStat,EStat,Tloc,Eloc)
	for j in range(9):
		T = TStat
		E = EStat
		if(j < 6):
			# print("inT")
			mecVar = np.zeros((6,1))
			if (Tvalue[i] == 0):
				delta = 100
			else:
				delta = -alph * Tvalue[i]
			mecVar[j] = delta

			T = TStat + mecVar
			S1,D1 = polyfunc(poly,T,E,Tinit + mecVar,Einit)[:2]

			T = TStat - mecVar
			S2,D2 = polyfunc(poly,T,E,Tinit - mecVar,Einit)[:2]

		else:
			elecVar = np.zeros((3,1))
			if(Evalue[k] == 0):
				delta = 100
			else:
				delta = alph * Evalue[k]
			elecVar[j-6] = delta

			E = EStat + elecVar
			S1,D1 = polyfunc(poly,T,E,Tinit,Einit + elecVar)[:2]

			E = EStat - elecVar
			S2,D2 = polyfunc(poly,T,E,Tinit,Einit - elecVar)[:2]

		# print(np.concatenate(((S1-S2)/(2*delta),(D1-D2)/(2*delta)),0).shape)
		toReturn.append((j,np.concatenate(((S1-S2)/(2*delta),(D1-D2)/(2*delta)),0)))
	return k,toReturn

main()