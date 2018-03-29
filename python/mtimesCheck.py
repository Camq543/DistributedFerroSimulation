import numpy as np
import numpy.matlib as npm
import time
import math
import mtimesx

workDir = '/home/lapis/Documents/FerroProject/python/DataFiles/mtimesCheck/'

def main():
	ranges = [10,50,math.pi,1e3,1.7e-7,9.43e5]
	steps = 6
	valueList = []
	vecList = []
	for val in ranges:
		vecList.append(np.linspace(0,val,steps))
	outer = []
	for i in range(len(vecList)):
		for j in range(len(vecList)):
			temp = npm.repmat(np.outer(vecList[i],vecList[j]),6,1).reshape(6,1,6,6)
			print(temp.shape)
			valueList.append(temp)
			

	for i in range(len(valueList)):
		data_write('n' + str(i + 1),mtimesx.mul(valueList[i],valueList[i],'N'))
		data_write('t' + str(i + 1),mtimesx.mul(valueList[i],valueList[i],'T'))
		data_write('f' + str(i + 1),mtimesx.mul(valueList[i],valueList[i],'F'))
	start = time.time()
	nVals = [10,50,100,1000,10000]
	mVals = [1,6,12,24]
	arrList = []
	outFile = open("/home/lapis/Documents/FerroProject/python/DataFiles/mtimesCheck/timespython.csv",'w')
	outFile.write('n,m,mode,time\n')
	for n in nVals:
		for m in mVals:
			temp = np.random.randn(n,m,6,6)
			times = np.zeros(10)
			for i in range(10):
				t1 = time.time()
				mtimesx.mul(temp,temp,'N')
				times[i] = time.time() - t1
			outFile.write(str(n) + ',' + str(m) + ',' + ',n,' + str(times.mean()) + '\n')
			times = np.zeros(10)
			for i in range(10):
				t1 = time.time()
				mtimesx.mul(temp,temp,'T')
				times[i] = time.time() - t1
			outFile.write(str(n) + ',' + str(m) + ',' + ',t,' + str(times.mean()) + '\n')
			times = np.zeros(10)
			for i in range(10):
				t1 = time.time()
				mtimesx.mul(temp,temp,'F')
				times[i] = time.time() - t1
			outFile.write(str(n) + ',' + str(m) + ',' + ',f,' + str(times.mean()) + '\n')
	outFile.close()

	

def data_write(fileName,toWrite):
	outFile = open(workDir + fileName + 'python.csv','w')
	outFile.write('grain,dir,i,j,data\n')
	for grain in range(len(toWrite)):
		for direction in range(len(toWrite[0])):
			for i in range(len(toWrite[0][0])):
				for j in range(len(toWrite[0][0][0])):
					outFile.write(str(grain) + ',' + str(direction) + ',' + str(i) + ',' + str(j) + ',' + str(toWrite[grain][direction][i][j]) + '\n')

main()
