#output processor
#R.M. Ward 2021

import json
import matplotlib.pyplot as plt
import numpy

SAVE_PLOTS = True

#load Output File
outputFileName = None
with open("CONFIG.txt", "r") as configFile:
    for line in configFile:
        if "ANALYSIS OUTPUT FILENAME" in line:
            outputFileName = line.split(":")[-1].strip()
            
if outputFileName is None:
    raise IOError("Could not locate filename from config.")

dataDict = None
with open(outputFileName, "r") as inFile:
    dataDict = json.load(inFile)
    
print("Data successfully loaded from "+str(outputFileName))
    
#first create data tables
for keyName in dataDict.keys():
    descriptor = "Full_Factorial" if keyName == "ff" else "Taguchi" if keyName == "tag" else "LHS" if keyName == "lhs" else "ERROR"
    caseList = dataDict[keyName]
    
    #h_L,fd1,fd2,oopThickness,massVal,maxMises
    outputText = descriptor+"\n"
    outputText += "t_skin [in],t_stiff [in],n_stiff [nondim.],h_stiff [in],w_stiff [in],mass [lb],critical load [lbf],sym/asym [nondim.]\n"
    for lineList in caseList:
        outputText += ",".join([str(x) for x in lineList])+"\n"
    
    outFileName = descriptor+"_outputTable.csv"
    with open(outFileName,"w") as outFile:
        outFile.write(outputText)
        
        
#write function to create factor effects
def plotFactorEffects(keyStr):
    try:
        descriptor = "Full_Factorial" if keyStr == "ff" else "Taguchi" if keyStr == "tag" else "ERROR"
        ffDat = dataDict[keyStr]
        numDV = 5
        numOV = 2
        
        dvNameList = ['t_skin','t_stiff','n_stiff','h_stiff','w_stiff']
        ovNameList = ['mass','crit_load']
        ovUnitList = ['[lb]','[lbf]']
        for dvIndex in range(numDV):
            #find possible DV values
            dvVals = []
            for i in range(len(ffDat)):
                dvVal = ffDat[i][dvIndex]
                if not dvVal in dvVals:
                    dvVals.append(dvVal)
            # print("Found values of "+dvNameList[dvIndex]+": "+str(dvVals))
            
            for ovIndex in range(numOV):
                ovIndexMapped = ovIndex+numDV
                #now find average OV of all values that include the dvVals
                foundOutputs = []
                for targetDVVal in dvVals:
                    listToAvg = []
                    for i in range(len(ffDat)):
                        if ffDat[i][dvIndex] == targetDVVal:
                            listToAvg.append(ffDat[i][ovIndexMapped])
                    # print(listToAvg)
                    avgVal = sum(listToAvg)/len(listToAvg)
                    foundOutputs.append(avgVal)
                    
                # print(foundOutputs)
                ylabelStr = ovNameList[ovIndex] + " " + ovUnitList[ovIndex]
                
                plt.xlim(0,len(dvVals)-1)
                plt.ylim(0,1000 if ovIndex == 0 else 2200000)
                plt.xlabel("Levels for DV= "+dvNameList[dvIndex])
                plt.ylabel(ylabelStr)
                plt.title(descriptor+" Factor Effects- DV: "+dvNameList[dvIndex]+"; OV: "+ovNameList[ovIndex])
                plt.plot(range(len(dvVals)),foundOutputs,marker=".")
                
                if SAVE_PLOTS:
                    plt.savefig(descriptor+"_DV-"+dvNameList[dvIndex]+"_OV-"+ovNameList[ovIndex]+".png")
                    plt.cla()
                else:
                    plt.show()
                
                
    except KeyError as kerr:
        print("Key not found: "+str(kerr))
    
# plotFactorEffects("ff")    
plotFactorEffects("tag")

#

