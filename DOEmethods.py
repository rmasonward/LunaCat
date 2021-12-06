"""
The Following script is to be used
for the structural optimization
course offered at Texas A&M University

It provides arrays non-dimensional values for 
each Design Variable based on given values 

Created by:
Adam Kellen '16
"""

import random
import numpy
import itertools

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def LHS(nV,nS):
    # Initialization:
    k=0
    # Creation of a list dictionary
    x = {}
    x[k] = []
    
    # Loop elements (part1)
    for i in range(nV):
    	x1 = []
    
    	for j in range(nS):
    		a = ((float(int(j)))/nS)
    		b = ((float(int(j))+1)/nS)
    		listesample = random.uniform(a,b)
    		x1.append(listesample)
    
    	# Select a random number nP times between each Sample and for each Var
    	for k in range(nS):
    		listechoice = random.choice(x1)
    		x.setdefault(k, []).append(listechoice)
    		x1.remove(listechoice)

    return list(x.values())

def FullFactorial(nV,nL):
             
    # Arrays for Itertools
    # temp 1 - possible answers    
    temp = [0.0]*nL
    for i in range(nL):
        temp[i] = i/(nL-1)

    # Itertools
    x = list(itertools.product(temp,repeat=nV))
    
    return x

def Taguchi(nV,FileName):
   
    DataFile = open(FileName)
    Data = numpy.loadtxt(DataFile)
    
    MaxVal = numpy.amax(Data)
    
    Fullx = (Data-1)/(MaxVal-1)
    
    x = Fullx[::, 0:nV]
    
    return x
