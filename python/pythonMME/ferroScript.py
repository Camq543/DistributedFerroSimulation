import numpy as np
import matplotlib.pyplot as plt
import time

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

	polycD = np.zeros((n,m))
	polycS = np.zeros((n,m))

	t1 = time.time()
	# outFile = open(workDir + "sim/" + "iterCount.csv",'w')
	# outFile.write("i,k,iterCount\n")
	for k in range(m):
		Tmacro = np.array([0,0,Tvalue[k],0,0,0]).reshape(6,1)

		Tloc = np.tile(Tmacro,(poly.nbg,1,1,1))
		Eloc = np.zeros((poly.nbg,1,3,1))

		for i in range(n):
			Emacro = np.array([0,0,Evalue[i]]).reshape(3,1)

			Smacro,Dmacro,Tloc,Eloc,iterCount = polyfunc(poly,Tmacro,Emacro,Tloc,Eloc)
			# outFile.write(str(i) + "," + str(k) + "," + str(iterCount) + "\n")

			polycD[i,k] = Dmacro[2,:]
			polycS[i,k] = Smacro[2,:]
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
			yvals.append(polycD[i,k])
		ax.plot(xvals,yvals,colors[k])

	fig2,ax2 = plt.subplots()
	for k in range(m):
		xvals = []
		yvals = []
		for i in range(n):
			xvals.append(Evalue[i])
			yvals.append(polycS[i,k] - polycS[0,k])
		ax2.plot(xvals,yvals,colors[k])

	#plt.show()

	nbt = Constants.nbt
	nbe = Constants.nbe
	alph = Constants.alpha
	LApp = np.zeros((nbt,nbe,9,9))
	Tvalue = np.linspace(0,Constants.tSpace,nbt) * 1e6
	#Tvalue = np.array([-75*1e6])
	Evalue = np.linspace(0,Constants.eSpace,nbe) * 1e6
	#Evalue = np.array([2.5*1e6])
	mecVar = np.zeros((6,1))
	elecVar = np.zeros((3,1))
	delta = 0

	t1 = time.time()
	for i in range(nbt):
		TStat = np.array([0,0,Tvalue[i],0,0,0]).reshape(6,1)
		for k in range(nbe):
			EStat = np.array([0,0,Evalue[k]]).reshape(3,1)
			# print(TStat)
			# print(EStat)
			Smacro,Dmacro,Tinit,Einit,iterCount = polyfunc(poly,TStat,EStat)
			for j in range(9):
				T = TStat
				E = EStat
				if(j < 6):
					mecVar = np.zeros((6,1))
					if (Tvalue[i] == 0):
						delta = 100
					else:
						delta = -alph * Tvalue[i]
					mecVar[j] = delta
					#print(mecVar)
					T = TStat + mecVar
					S1,D1 = polyfunc(poly,T,E,Tinit + np.tile(mecVar,(poly.nbg,1,1,1)),Einit)[:2]

					T = TStat - mecVar
					S2,D2 = polyfunc(poly,T,E,Tinit - np.tile(mecVar,(poly.nbg,1,1,1)),Einit)[:2]

				else:
					elecVar = np.zeros((3,1))
					if(Evalue[k] == 0):
						delta = 100
					else:
						delta = alph * Evalue[k]
					elecVar[j-6] = delta
					#print(elecVar)
					E = EStat + elecVar
					S1,D1 = polyfunc(poly,T,E,Tinit,Einit + np.tile(elecVar,(poly.nbg,1,1,1)))[:2]

					E = EStat - elecVar
					S2,D2 = polyfunc(poly,T,E,Tinit,Einit - np.tile(elecVar,(poly.nbg,1,1,1)))[:2]

				temp = np.concatenate(((S1-S2)/(2*delta),(D1-D2)/(2*delta)),0)

				LApp[i,k,:,j] = temp[:,0]

	write_csv("LApp",LApp)

	piezoTime = time.time() - t1
	print("Piezo coefficients took " + str(piezoTime) + " seconds.")



main()