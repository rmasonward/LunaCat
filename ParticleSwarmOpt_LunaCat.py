#Particle Swarm Optimizer
#R. Mason Ward

import random
import DOEmethods as doe
from LunaCatAbqMainCode import evalModel


# random.seed(50)
#Pick final result output file to append
optListFileName = "OPT_LOG.csv"
runListFileName = "RUNNING_LOG.csv"
currentBestFileName = "CURRENT_BEST.txt"

#Configuration
POP_SIZE = 6
MAX_ITERS = 10

NUM_VARS = 11

arm_base_width_min = 0.1 #m #
arm_base_width_max = 0.4 #m

arm_base_height_min = 0.1 #m #
arm_base_height_max = 0.3 #m

arm_wall_thickness_min = 0.001 #m #
arm_wall_thickness_max = 0.0225 #m

arm_length_min = 2.0 #m
arm_length_max = 5.0 #m #

arm_taper_ratio_min = 0.5 #
arm_taper_ratio_max = 1.0

wall_length_min = 2.0 #m
wall_length_max = 4.0 #m #

axle_length_min = 2.0 #m #
axle_length_max = 4.0 #m

#KDVs
top_opt_thickness_min = 0.05 #m #
top_opt_thickness_max = 0.2 #m

L1_percent_min = 5
L1_percent_max = 20

# L2_percent_min = 5
# L2_percent_max = 20

clevis_edge_thickness_min = 0.01 #
clevis_edge_thickness_max = 0.1

# side_wall_thickness_min = 0.01 #
# side_wall_thickness_max = 0.5

#Axle/arm material (discrete DV)
material_type_min = 1
material_type_max = 3 #

######## DEFINE YEILD STRESSES #########
mat1_yield = 324e6
mat2_yield = 276e6
mat3_yield = 827e6

######### DEFINE OBJECTIVES #########
V_GOAL = 70.55 #m/s
ANGLE_GOAL = 50.81 #deg
MASS_REF = 5029.582494 #kg #This value is used to non-dimensionalize the mass in the cost function; currently chosen to the structure set to the upper bound of all DVs: 0.4, 0.3, 0.0225, 5.0, 1.0, 4.0, 4.0, 3, 0.2, 20.0, 0.1

#t_skin [m],t_stiff [m],n_stiff [nondim.],h_stiff [m],w_stiff [m]
boundsList = [[arm_base_width_min, arm_base_width_max], [arm_base_height_min, arm_base_height_max], [arm_wall_thickness_min, arm_wall_thickness_max], [arm_length_min, arm_length_max], 
[arm_taper_ratio_min, arm_taper_ratio_max], [wall_length_min, wall_length_max], [axle_length_min, axle_length_max], [material_type_min, material_type_max], [top_opt_thickness_min, top_opt_thickness_max], 
[L1_percent_min, L1_percent_max],[clevis_edge_thickness_min, clevis_edge_thickness_max]]

#Tunable Params
V_MAX = 0.5
ALPHA = 0.25#0.3
P_BEST_COEFF = 0.85
G_BEST_COEFF = 0.025#0.2

numEvals = 1

#To be used for convergence plotting later
iterNums = []
massNums = []

#[t_skin,t_stiff,n_stiff,h_stiff,w_stiff,gFuncBest]
with open(currentBestFileName, "w") as outFile:
    outFile.write("Global Best Update Log\n")

with open(optListFileName, "w") as outFile:
    outFile.write("OPTIMIZATION LOG-BEST PER GENERATION\narm_base_width,arm_base_height,arm_wall_thickness,arm_length,arm_taper_ratio,wall_length,axle_length,material_type,top_opt_thickness,L1,clevis_edge_thickness,costVal\n")
    
with open(runListFileName, "w") as outFile:
    outFile.write("RUN LOG\narm_base_width,arm_base_height,arm_wall_thickness,arm_length,arm_taper_ratio,wall_length,axle_length,material_type,top_opt_thickness,L1,clevis_edge_thickness,veloMagMax,veloAngle,eigenVal1,maxMises,mass,costVal\n")    

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
    arm_base_width,arm_base_height,arm_wall_thickness,arm_length,arm_taper_ratio,wall_length,axle_length,material_type,top_opt_thickness,L1,clevis_edge_thickness = mapNormedVecToFullScale(varVec)
    material_type = int(round(material_type)) #discretize material type
    
    global numEvals
    
    #Now compute cost func and buckling criteria
    try:
        veloMagMax,veloAngle,eigenVal1,maxMises,mass = evalModel(arm_base_width,arm_base_height,arm_wall_thickness,arm_length,arm_taper_ratio,wall_length,axle_length,material_type,top_opt_thickness,L1,clevis_edge_thickness,numEvals) #Run Abq Functional Evaluation
        
        doesBuckle = eigenVal1 <= 0 #buckling constraint from problem
        doesExceedMises = maxMises >= (mat1_yield if material_type == 1 else mat2_yield if material_type == 2 else mat3_yield)
        
        #violates constraints?
        doesViolateConstraints = doesBuckle or doesExceedMises
        
        #convert angles back to radian to nondimensionalize
        veloAngleRad = veloAngle*np.pi/180
        angleGoalRad = ANGLE_GOAL*np.pi/180
        
        #compute cost value
        costVal = 2*(veloMagMax-V_GOAL)**2/V_GOAL**2 + (veloAngle-ANGLE_GOAL)**2 + mass/MASS_REF
        
        #apply penalty constraint
        costVal = 1e15 if doesViolateConstraints else costVal
    except Exception as e:
            print("EXCEPTION: "+str(e))
            veloMagMax=-1e20
            veloAngle=0.0
            eigenVal1=-1e20
            maxMises=1e20
            mass=1e20
            costVal = 1e20
        
    numEvals += 1
    
    
    #Add best after each generation
    
    fullScaleOutputList = [arm_base_width,arm_base_height,arm_wall_thickness,arm_length,arm_taper_ratio,wall_length,axle_length,material_type,top_opt_thickness,L1,clevis_edge_thickness,veloMagMax,veloAngle,eigenVal1,maxMises,mass,costVal]
    print("DV vector full-scale: "+str(fullScaleOutputList))
    
    with open(runListFileName, "a") as outFile:
        outFile.write(",".join([str(x) for x in fullScaleOutputList]) +"\n")
        
        
    return costVal
    

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
                
                #Rewrite file if new glob best
                arm_base_width,arm_base_height,arm_wall_thickness,arm_length,arm_taper_ratio,wall_length,axle_length,material_type,top_opt_thickness,L1,clevis_edge_thickness = mapNormedVecToFullScale(gXBest)
                material_type = int(round(material_type))
                fullScaleOutputList = [arm_base_width,arm_base_height,arm_wall_thickness,arm_length,arm_taper_ratio,wall_length,axle_length,material_type,top_opt_thickness,L1,clevis_edge_thickness,gFuncBest]
                print("DV vector full-scale: "+str(fullScaleOutputList))

                with open(currentBestFileName, "a") as outFile:
                    outFile.write(",".join([str(x) for x in fullScaleOutputList]) +",Job-"+str(numEvals-1)+"\n")
                
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
                
                #Rewrite file if new glob best
                arm_base_width,arm_base_height,arm_wall_thickness,arm_length,arm_taper_ratio,wall_length,axle_length,material_type,top_opt_thickness,L1,clevis_edge_thickness = mapNormedVecToFullScale(gXBest)
                material_type = int(round(material_type))
                fullScaleOutputList = [arm_base_width,arm_base_height,arm_wall_thickness,arm_length,arm_taper_ratio,wall_length,axle_length,material_type,top_opt_thickness,L1,clevis_edge_thickness,gFuncBest]
                print("DV vector full-scale: "+str(fullScaleOutputList))

                with open(currentBestFileName, "a") as outFile:
                    outFile.write(",".join([str(x) for x in fullScaleOutputList]) +",Job-"+str(numEvals-1)+"\n")
        
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
    arm_base_width,arm_base_height,arm_wall_thickness,arm_length,arm_taper_ratio,wall_length,axle_length,material_type,top_opt_thickness,L1,clevis_edge_thickness = mapNormedVecToFullScale(gXBest)
    material_type = int(round(material_type))
    fullScaleOutputList = [arm_base_width,arm_base_height,arm_wall_thickness,arm_length,arm_taper_ratio,wall_length,axle_length,material_type,top_opt_thickness,L1,clevis_edge_thickness,gFuncBest]
    print("DV vector full-scale: "+str(fullScaleOutputList))

    with open(optListFileName, "a") as outFile:
        outFile.write(",".join([str(x) for x in fullScaleOutputList]) +"\n")

#now extract global best after specified number of iterations
# print(costFunc(gXBest),gFuncBest)
print("Minimum cost: "+str(gFuncBest))
print("DV vector: "+str(gXBest))

print("DV vector full-scale: "+str(mapNormedVecToFullScale(gXBest)))

print("Done.")

