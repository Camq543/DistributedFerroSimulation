#python/3rd party libs
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing
import time
import os
from functools import partial
from pyspark import SparkContext
from pyspark import SparkConf

#my functions
from runFunctions import run_poly, piezo_map_nbt, mongo_data_fill
from util_functions.matToFile import write_csv
from constants import Constants
from polyfunc import polyfunc
from structures.monocrystal import MonoCrystal as Monoc
from structures.polycrystal import PolyCrystal as Polyc

def main():
	try:
		startTime = time.time()
		workDir = "/home/lapis/Documents/FerroProject/python/pythonMME/DataFiles/"
		mono = Monoc()
	#	mono.data_write()
		poly = Polyc(mono)
	#	poly.data_write()

		polyTime = run_poly(poly)
		print("Polycrystal function took " + str(polyTime) + " seconds.")
		#plot_polycrystal(polycS,polycD,Evalue,Tvalue)

		sconf = SparkConf().set("spark.default.parallelism",int(os.environ['SLURM_JOB_NUM_NODES']))
		sconf = sconf.set("spark.local.dir","$HOME/tmp/spark-logs/")
		sconf = sconf.set("spark.locality.wait",10)
		master_node = 'spark://' + os.environ['MASTER'] + ':7077'
		
		
		sc = SparkContext(master = master_node,conf = sconf)
		time.sleep(5)
		t1 = time.time()
		#poly = sc.broadcast(poly)
		LApp = piezo_map_nbt(poly,sc)
		piezoTime = time.time() - t1

		#write_csv("LApp",LApp)

		totalTime = time.time()-startTime
		print("Piezo coefficents took " + str(piezoTime) + " seconds.")
		mongo_data_fill(polyTime,piezoTime,totalTime)
		time.sleep(5)
		print('Cancelling')
		os.system("scancel %s" % os.environ['SLURM_JOBID'])
	except Exception as e:
		print("EXCEPTINO THROWN, CANCELLING")
		os.system("scancel %s" % os.environ['SLURM_JOBID'])
		raise e

def plot_polycrystal(polycS, polycD, Evalue, Tvalue):
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

	plt.show()

main()