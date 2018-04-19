#write tensor to a csv file
def write_csv(filename,toWrite):
	outFile = open(filename + 'python.csv','w')
	shape = toWrite.shape
	if len(shape) == 4:
		outFile.write('grain,dir,i,j,data\n')
		for grain in range(shape[0]):
			for direc in range(shape[1]):
				for i in range(shape[2]):
					for j in range(shape[3]):
						outFile.write(','.join((str(grain),str(direc),str(i),str(j),str(toWrite[grain,direc,i,j]))) + '\n')
		outFile.close()
	elif len(shape) == 5:
		outFile.write('grain,dir,i,j,k,data\n')
		for grain in range(shape[0]):
			for direc in range(shape[1]):
				for i in range(shape[2]):
					for j in range(shape[3]):
						for k in range(shape[4]):
							outFile.write(','.join((str(grain),str(direc),str(i),str(j),str(k),str(toWrite[grain,direc,i,j,k]))) + '\n')
		outFile.close()
	elif len(shape) == 6:
		outFile.write('grain,dir,i,j,k,l,data\n')
		for grain in range(shape[0]):
			for direc in range(shape[1]):
				for i in range(shape[2]):
					for j in range(shape[3]):
						for k in range(shape[4]):
							for l in range(shape[5]):
								outFile.write(','.join((str(grain),str(direc),str(i),str(j),str(k),str(l),str(toWrite[grain,direc,i,j,k,l]))) + '\n')
		outFile.close()		
