#Particle Swarm Optimizer
#R. Mason Ward

import random
import DOEmethods as doe
from LunaCatAbqMainCode import evalModel


# random.seed(50)
#Pick final result output file to append
optListFileName = "OPT_LOG.csv"

#Configuration
POP_SIZE = 14
NUM_VARS = 5
MAX_ITERS = 12

t_skin_min = 0.025
t_skin_max = 0.25

t_stiff_min = 0.025
t_stiff_max = 0.25

n_stiff_min = 3
n_stiff_max = 10

h_stiff_min = 0.25
h_stiff_max = 5.0 

w_stiff = 0.15

SEED_SIZE = 0.5

#width of the stiffener was removed with DOE 
# w_stiff_min = 0.25
# w_stiff_max = 0.25 #chosen so that for 10 stiffeners per 60 inches, the flanges do not overlap

#t_skin [m],t_stiff [m],n_stiff [nondim.],h_stiff [m],w_stiff [m]
boundsList = [[t_skin_min,t_skin_max], [t_stiff_min,t_stiff_max], [n_stiff_min,n_stiff_max], [h_stiff_min,h_stiff_max]]

#Tunable Params
V_MAX = 0.5
ALPHA = 0.25#0.3
P_BEST_COEFF = 0.85
G_BEST_COEFF = 0.025#0.2

#To be used for convergence plotting later
iterNums = []
massNums = []

#[t_skin,t_stiff,n_stiff,h_stiff,w_stiff,gFuncBest]
with open(optListFileName, "w") as outFile:
    outFile.write(",".join(['t_skin','t_stiff','n_stiff','h_stiff','w_stiff','mass']) +"\n")

def mapNormedVecToFullScale(varVec):
    listToAdd = []
    for i,boundTuple in enumerate(boundsList):
        percentVal = varVec[i]
        mappedVal = (boundTuple[1]-boundTuple[0])*percentVal + boundTuple[0]
        listToAdd.append(mappedVal)
    
    return listToAdd

#objective function
def costFunc(varVec):
    #this function should extract specific objective values and return a single cost value to be minimized
    
    #compute an unnormalized DV vector to be supplied to the cost function evaluation
    t_skin,t_stiff,n_stiff,h_stiff = mapNormedVecToFullScale(varVec)
        
    n_stiff = int(round(n_stiff)) #Bucket number of stiffeners into discrete values
    
    #Now compute cost func and buckling criteria
    massVal,critForce,symStr = evalModel(n_stiff, t_stiff, h_stiff, w_stiff, t_skin, SEED_SIZE)
    
    w_domain = 60./n_stiff
    doesBuckle = critForce/w_domain <= 40 #buckling constraint from problem
    
    #apply penalty constraint
    massVal = 1e15 if doesBuckle else massVal
    # if massVal == 1e15:
        # print("BUCKLES!")
        # print(t_skin,t_stiff,n_stiff,h_stiff,w_stiff)
    # else:
        # print("NO BUCKLE! -> "+str(massVal))
        # print(varVec)
    return massVal
    

#Data structure initialization
populationList = []

gFuncBest = 1e15
gXBest = [None for i in range(NUM_VARS)]

indivFuncBest = [1e15 for indiv in range(POP_SIZE)]
indivXBest = [[None for i in range(NUM_VARS)] for indiv in range(POP_SIZE)] 
indivVelos = [[0 for i in range(NUM_VARS)] for indiv in range(POP_SIZE)]

#seed initial variable vector    
#LHS implemented here.
populationList = doe.LHS(NUM_VARS,POP_SIZE)
    
#now begin iterative loop
iterCnt = 0
while iterCnt <= MAX_ITERS:
    iterCnt += 1
    print("####################\nBeginning Iteration "+str(iterCnt)+"/"+str(MAX_ITERS)+"\n####################\n")
    
    if iterCnt == 1: #seed initial G Best
        for indivIndex in range(POP_SIZE):
            #extract pertinent individual values
            indivXVector = populationList[indivIndex][:]
            
            indivFuncVal = costFunc(indivXVector)

            indivFuncBest[indivIndex] = indivFuncVal
            indivXBest[indivIndex] = indivXVector
                
            #now check global best too
            if indivFuncVal < gFuncBest:
                gFuncBest = indivFuncVal
                gXBest = indivXVector[:]
                
        # print("ITER 1: "+str(gFuncBest)+"; "+str(gXBest))
    
    nextGFuncBest = gFuncBest
    nextGXBest = gXBest[:]
        
    for indivIndex in range(POP_SIZE):
        #extract pertinent individual values
        indivXVector = populationList[indivIndex][:]
        
        indivFuncVal = costFunc(indivXVector)
        
        #compare to personal best
        if indivFuncVal < indivFuncBest[indivIndex]:
            indivFuncBest[indivIndex] = indivFuncVal
            indivXBest[indivIndex] = indivXVector
            
            #now check global best too
            if indivFuncVal < nextGFuncBest:
                nextGFuncBest = indivFuncVal
                nextGXBest = indivXVector[:]
                # print(indivFuncVal,costFunc(indivXVector))
        
        #compute velocity to new point
        for dimIndex in range(NUM_VARS):
            lowBound = 0
            upBound = 1 #boundsList[dimIndex]
        
            newComponentVelo = ALPHA*indivVelos[indivIndex][dimIndex] + P_BEST_COEFF*random.random()*(indivXBest[indivIndex][dimIndex]-indivXVector[dimIndex]) + G_BEST_COEFF**random.random()*(gXBest[dimIndex]-indivXVector[dimIndex])
            newComponentVelo = min(newComponentVelo,V_MAX) # restrict velocity component to maximum if it exceeds
            
            newXComponent = indivXVector[dimIndex] + newComponentVelo
            newXComponent = min(max(newXComponent,lowBound),upBound) #restrict to bound
            
            #update array values
            populationList[indivIndex][dimIndex] = newXComponent
            indivVelos[indivIndex][dimIndex] = newComponentVelo
        
    gFuncBest = nextGFuncBest
    gXBest = nextGXBest[:]
    
    iterNums.append(iterCnt)
    massNums.append(gFuncBest)
    
    #Add best after each generation
    t_skin,t_stiff,n_stiff,h_stiff = mapNormedVecToFullScale(gXBest)
    n_stiff = int(round(n_stiff))
    fullScaleOutputList = [t_skin,t_stiff,n_stiff,h_stiff,w_stiff,gFuncBest]
    print("DV vector full-scale: "+str(fullScaleOutputList))

    with open(optListFileName, "a") as outFile:
        outFile.write(",".join([str(x) for x in fullScaleOutputList]) +"\n")

#now extract global best after specified number of iterations
# print(costFunc(gXBest),gFuncBest)
print("Minimum cost: "+str(gFuncBest))
print("DV vector: "+str(gXBest))

print("DV vector full-scale: "+str(mapNormedVecToFullScale(gXBest)+[w_stiff]))

print("Done.")

