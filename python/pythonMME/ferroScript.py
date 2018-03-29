import numpy as np
import matplotlib.pyplot as plt
from polyfunc import polyfunc
from structures.monocrystal import MonoCrystal as Monoc
from structures.polycrystal import PolyCrystal as Polyc

def main():
	mono = Monoc()
	mono.data_write()
	poly = Polyc(mono)
	poly.data_write()

	n = 51
	#n = 1
	Evalue = np.linspace(0,2e6,n)
	Tvalue = np.array([-100,-50,0,50,100]) * 1e6
	#Evalue = np.array([1e6])
	#Tvalue = np.array([10e7])
	m = len(Tvalue)

	polycD = np.zeros((n,m))
	polycS = np.zeros((n,m))

	for k in range(m):
		Tmacro = np.array([0,0,Tvalue[k],0,0,0]).reshape(6,1)

		Tloc = np.tile(Tmacro,(poly.nbg,1,1,1))
		Eloc = np.zeros((poly.nbg,1,3,1))

		for i in range(n):
			#print("i: " + str(i) + "k: " + str(k))
			Emacro = np.array([0,0,Evalue[i]]).reshape(3,1)
			#print(Tmacro)
			#print(Emacro)
			Smacro,Dmacro,Tloc,Eloc,iterCount = polyfunc(poly,Tmacro,Emacro,Tloc,Eloc)

			polycD[i,k] = Dmacro[2,:]
			polycS[i,k] = Smacro[2,:]

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


	plt.show()





main()