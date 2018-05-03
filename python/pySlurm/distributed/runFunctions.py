#python/3rd party libs
import os
from os.path import expanduser
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing
import time
from pymongo import MongoClient
from functools import partial
from pyspark import SparkContext

from util_functions.matToFile import write_csv
from constants import Constants
from polyfunc import polyfunc
from structures.monocrystal import MonoCrystal as Monoc
from structures.polycrystal import PolyCrystal as Polyc

def run_poly(poly):
    ##Runs initial round of iterations to determine the behavior of the polycrystal
    n = 51
    Evalue = np.linspace(0,2e6,n) #values to use for electric field
    Tvalue = np.array([-100,-50,0,50,100]) * 1e6 #values to use for mechanical stress
    m = len(Tvalue)

    t1 = time.time() #start timer
    pool = multiprocessing.Pool() #prepare parallel pool
    func = partial(map_poly,poly,Tvalue,Evalue) #prepare function for parallelization
    results = pool.map(func, range(m)) #get results from parallel computation
    pool.close() #close parallel pool

    #initialize arrays to hold polycrystal behavior
    polycD = np.zeros((m,n))
    polycS = np.zeros((m,n))

    #store values obtained in polycrystal calculation for use plotting
    for k in range(m):
        polycD[k,:] = results[k][1][:]
        polycS[k,:] = results[k][2][:]

    return time.time() - t1, polycS, polycD #return computation time

def map_poly(poly,Tvalue,Evalue,k):
    #initialize arrays to hold results
    tempD = np.zeros(len(Evalue))
    tempS = np.zeros(len(Evalue))

    Tmacro = np.array([0,0,Tvalue[k],0,0,0]).reshape(6,1) #create stress vector based on input

    #create tensors to be used in polyfunc
    Tloc = np.tile(Tmacro,(poly.nbg,1,1,1))
    Eloc = np.zeros((poly.nbg,1,3,1))

    #iterate over different Evalues and determine polycrystal behavior
    for i in range(len(Evalue)):
        Emacro = np.array([0,0,Evalue[i]]).reshape(3,1)

        Smacro,Dmacro,Tloc,Eloc,iterCount = polyfunc(poly,Tmacro,Emacro,Tloc,Eloc)

        tempD[i] = Dmacro[2,:]
        tempS[i] = Smacro[2,:]

    return k,tempD,tempS #return k value given along with stress and electric behavior

def piezo_map_tasks(poly,sc):
    nbt = Constants.nbt
    nbe = Constants.nbe
    LApp = np.zeros((nbt,nbe,9,9)) #array to store results
    times = np.zeros((nbt,nbe)) #array to store computation times
    taskArr = [] #array to partition into tasks

    poly = sc.broadcast(poly)
    Tvalue = sc.broadcast(np.linspace(0,Constants.tSpace,nbt) * 1e6) #vector of stress values to use
    Evalue = sc.broadcast(np.linspace(0,Constants.eSpace,nbe) * 1e6) #vector of electric values to use
    #add all index pairs used to pull values from Evalue and Tvalue
    for i in range(nbt):
        for k in range(nbe):
            taskArr.append((i,k))
    
    print('DEFAULT PARALLELISM',sc.defaultParallelism)
    tMap = sc.parallelize(taskArr,nbt*nbe) #parallelize across all index pairs
    print("NUMBER OF PARTITIONS", tMap.getNumPartitions())
    func = partial(run_piezo_task,poly,Tvalue,Evalue) #prepare function for parallelization
    tMap = tMap.map(lambda i: func(i)) #create tasks for each index pair in taskArr
    tMap = tMap.collect() #collect results
    
    #unpack data and store results
    for tDat in tMap:
        i,k = tDat[0] 
        times[i,k] = tDat[2] #store computation time
        for jDat in tDat[1]:
            LApp[i,k,:,jDat[0]] = jDat[1][:,0] #store results
    
    return LApp, times #return results and computation times

def piezo_map_nbt(poly,sc):
    #map stress values to be used for computation
    nbt = Constants.nbt
    nbe = Constants.nbe
    LApp = np.zeros((nbt,nbe,9,9)) #array to store results

    print('DEFAULT PARALLELISM',sc.defaultParallelism) 
    tMap = sc.parallelize(range(nbt),nbt) #distribute across a range of stress values represented by indexes in range(nbt)
    print("NUMBER OF PARTITIONS", tMap.getNumPartitions())
    tMap = tMap.map(lambda i: map_piezo_e(poly,i)) #create task using poly for each stress index
    tMap = tMap.collect() #collect function results
    
    #store data in LApp array
    for tDat in tMap:
        for eDat in tDat[1]:
            for jDat in eDat[1]:
                LApp[tDat[0],eDat[0],:,jDat[0]] = jDat[1][:,0]
    
    return LApp

def map_piezo_e(poly,i):
    #for a given stress index, parallelize across different electric fields
    nbe = Constants.nbe
    nbt = Constants.nbt
    Tvalue = np.linspace(0,Constants.tSpace,nbt) * 1e6 #vector of stress values to use
    Evalue = np.linspace(0,Constants.eSpace,nbe) * 1e6 #vector of electric field values to use

    pool = multiprocessing.Pool(processes = int(os.environ['SLURM_CPUS_PER_TASK'])) #initialize parallel pool based on number of cpus needed per task
    func = partial(run_piezo, poly, Tvalue, Evalue,i) #prepare function for parallelization
    data = pool.map(func,range(nbe))#parallelize across a range electric field values represented by indexes in range(nbe)
    pool.close()
    return i,data #return stress index and calculated data

def run_piezo(poly,Tvalue,Evalue,i,k):
    #determine piezo coefficients for given polycrystal 
    alph = Constants.alpha #determines degree of variation in polyfunc tests
    delta = 0 #initialize var to hold degree of variation 
    TStat = np.array([0,0,Tvalue[i],0,0,0]).reshape(6,1) #stress vector for use in polyfunc
    EStat = np.array([0,0,Evalue[k]]).reshape(3,1) #electric field vector for use in polyfunc
    Tloc = np.tile(TStat,(poly.nbg,1,1,1)) #stress tensor for use in polyfunc
    Eloc = np.tile(EStat,(poly.nbg,1,1,1)) #electric field tensor for use in polyfunc
    toReturn = []

    Smacro,Dmacro,Tinit,Einit,iterCount = polyfunc(poly,TStat,EStat,Tloc,Eloc) #calculate starting tensors for use below
    
    #stress vectors are 6x1, electric field are 3x1
    #we iterate over the 6 possible stress indexes, and 3 possible electric indexes, varying one each time
    #In the end, we have a 9x9 vector created from the stress and induction response of the polycrystal based on the 9 different variations
    for j in range(9): 
        T = TStat
        E = EStat
        if(j < 6):
            #if j <6, we vary stress
            mecVar = np.zeros((6,1)) #initialize variance array
            if (Tvalue[i] == 0):
                delta = 100
            else:
                delta = -alph * Tvalue[i]
            mecVar[j] = delta #set index j to a value based on alpha and our given stress values 
            #add variation to original stress vector and calculate
            T = TStat + mecVar 
            S1,D1 = polyfunc(poly,T,E,Tinit + mecVar,Einit)[:2]
            #subtract variation from original electric field vector and calculate
            T = TStat - mecVar
            S2,D2 = polyfunc(poly,T,E,Tinit - mecVar,Einit)[:2]

        else:
            #if j > 6, we vary electric field, logic is identical to above calculations
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

        toReturn.append((j,np.concatenate(((S1-S2)/(2*delta),(D1-D2)/(2*delta)),0))) #concatenate stress and induction response, and divide by 2*delta
    return k,toReturn #return given electric field index and response array

def run_piezo_task(poly,Tvalue,Evalue,indexes):
    #determine piezo coefficients for given polycrystal 
    #as opposed to the run_piezo function, index pairs are passed in as a tuple to allow dynamic allocation of tasks
    
    t1 = time.time() 
    alph = Constants.alpha #determines degree of variation in polyfunc tests
    delta = 0 #initialize var to hold degree of variation 
    i,k = indexes #unpack indexes
    TStat = np.array([0,0,Tvalue.value[i],0,0,0]).reshape(6,1) #stress vector for use in polyfunc
    EStat = np.array([0,0,Evalue.value[k]]).reshape(3,1) #electric field vector for use in polyfunc
    Tloc = np.tile(TStat,(poly.value.nbg,1,1,1)) #stress tensor for use in polyfunc
    Eloc = np.tile(EStat,(poly.value.nbg,1,1,1)) #electric field tensor for use in polyfunc
    toReturn = []

    Smacro,Dmacro,Tinit,Einit,iterCount = polyfunc(poly.value,TStat,EStat,Tloc,Eloc) #calculate starting tensors for use below
    
    #stress vectors are 6x1, electric field are 3x1
    #we iterate over the 6 possible stress indexes, and 3 possible electric indexes, varying one each time
    #In the end, we have a 9x9 vector created from the stress and induction response of the polycrystal based on the 9 different variations
    for j in range(9):
        T = TStat
        E = EStat 
        if(j < 6):
            #if j <6, we vary stress
            mecVar = np.zeros((6,1)) #initialize variance array
            if (Tvalue.value[i] == 0):
                delta = 100
            else:
                delta = -alph * Tvalue.value[i]
            mecVar[j] = delta #set index j to a value based on alpha and our given stress values 
            #add variation to original stress vector and calculate
            T = TStat + mecVar 
            S1,D1 = polyfunc(poly.value,T,E,Tinit + mecVar,Einit)[:2]
            #subtract variation from original electric field vector and calculate
            T = TStat - mecVar
            S2,D2 = polyfunc(poly.value,T,E,Tinit - mecVar,Einit)[:2]

        else:
            #if j > 6, we vary electric field, logic is identical to above calculations
            elecVar = np.zeros((3,1))
            if(Evalue.value[k] == 0):
                delta = 100
            else:
                delta = alph * Evalue.value[k]
            elecVar[j-6] = delta

            E = EStat + elecVar
            S1,D1 = polyfunc(poly.value,T,E,Tinit,Einit + elecVar)[:2]

            E = EStat - elecVar
            S2,D2 = polyfunc(poly.value,T,E,Tinit,Einit - elecVar)[:2]

        toReturn.append((j,np.concatenate(((S1-S2)/(2*delta),(D1-D2)/(2*delta)),0))) #concatenate stress and induction response, and divide by 2*delta
    return indexes,toReturn,time.time() - t1 #return given indexes, response array, and computation time

def mongo_data_fill(polycTime, piezoTime, totalTime,taskTimes):
    #add calculation time to mongo database
    #----------------------MONGODB STUFF------------------------------------#
    with open(expanduser('~/.ipmongo')) as f:
        IPadd = f.readline().strip()
    client = MongoClient(IPadd, 27017)

    db = client['ferro-runs']
    pycollection = db['pycollection']
    breakcollection = db['pybreak']
    
    grains = Constants.n1*Constants.n2*Constants.n3
    insert_data = {'num-nodes':int(os.environ['SLURM_JOB_NUM_NODES']), 'polyc-time':polycTime, 'piezo-time':piezoTime, 'total-time':totalTime,'nbt':Constants.nbt,'nbe':Constants.nbe,'grains':grains}
    if Constants.tobreak:
        breakcollection.insert_one(insert_data)
    else:
        pycollection.insert_one(insert_data)

    taskCollection = db['tasktimes']
    taskbreak = db['taskbreak']
    for i in range(Constants.nbt):
        for j in range(Constants.nbe):
            insert_data = {'i':i,'k':j,"nbt":Constants.nbt,"nbe":Constants.nbe,"time":taskTimes[i,j],'grains':grains,'num-nodes':int(os.environ['SLURM_JOB_NUM_NODES'])}
            if Constants.tobreak:
                taskbreak.insert_one(insert_data)
            else:
                taskCollection.insert_one(insert_data)