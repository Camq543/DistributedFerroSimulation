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
from runFunctions import run_poly, piezo_map_nbt, mongo_data_fill, piezo_map_tasks
#from util_functions.matToFile import write_csv
from constants import Constants
from polyfunc import polyfunc
from structures.monocrystal import MonoCrystal as Monoc
from structures.polycrystal import PolyCrystal as Polyc

def main():
	try:
		startTime = time.time()
		workDir = "/home/lapis/Documents/FerroProject/python/pythonMME/DataFiles/"
		mono = Monoc() #initialize monocrystal
	#	mono.data_write()
		poly = Polyc(mono) #using the monocrystal, initialize polycrystal
	#	poly.data_write()

		polyTime, polycS, polycD = run_poly(poly) #determine polycrystal behavior
		print("Polycrystal function took " + str(polyTime) + " seconds.")
		#plot_polycrystal(polycS,polycD,Evalue,Tvalue) #plot polycrystal behavior

		#set up configuration for spark, namely logging directory and default parallelism
		#we want to treat each piezo calculation as one task, so we parallelize across each combination of range(nbt) and range(nbe)
		sconf = SparkConf().set("spark.default.parallelism",Constants.nbt * Constants.nbe)
		#sconf = sconf.set("spark.local.dir","$HOME/tmp/spark-logs/")
		#master_node = 'spark://' + os.environ['MASTER'] + ':7077' #keep track of where the master node is
		
		#initialize SparkContext, which is used to interface between python and spark
		sc = SparkContext(master = "local[4]",conf = sconf) #when using Slurm, we submit with master = master_node, but on local pc just use local
		time.sleep(5)#wait until spark is finished setting up 
		t1 = time.time()
		LApp,taskTimes = piezo_map_tasks(poly,sc) #pass our polycrystal and SparkContext for use in calculation of piezo coefficients
		piezoTime = time.time() - t1

		totalTime = time.time()-startTime
		print("Piezo coefficents took " + str(piezoTime) + " seconds.")
		#mongo_data_fill(polyTime,piezoTime,totalTime,taskTimes) #store our run data in mongo db
		#print('Cancelling')
		#os.system("scancel %s" % os.environ['SLURM_JOBID']) #cancel slurm job to shutdown cluster
	except Exception as e:
		print("EXCEPTION THROWN, CANCELLING")
		#print(e)
		#os.system("scancel %s" % os.environ['SLURM_JOBID']) #in case of error, cancel slurm job
		raise e

def plot_polycrystal(polycS, polycD, Evalue, Tvalue):
	#plot behavior of polycrystal
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