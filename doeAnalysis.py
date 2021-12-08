#Main analysis integration code for DOE

import json
from LunaCatAbqMainCode import evalModel


inputFileName = None
outputFileName = None
with open("CONFIG.txt", "r") as configFile:
    for line in configFile:
        if "INPUT FILENAME" in line:
            inputFileName = line.split(":")[-1].strip()
        elif "ANALYSIS OUTPUT FILENAME" in line:
            outputFileName = line.split(":")[-1].strip()
            
if inputFileName is None or outputFileName is None:
    raise IOError("Could not locate filenames from config.")
    
#Load data from JSON
dataDict = None
with open(inputFileName, "r") as inFile:
    dataDict = json.load(inFile)
    
print("Data successfully loaded from "+str(inputFileName))

#loop thru data and perform functional evaluations
outputDict = {}
with open("goodBad.txt", "w") as outFile:
    outFile.write("LOG\n________________________\n")
    
with open("safetyNet.txt", "w") as outFile:
    outFile.write("arm_base_width,arm_base_height,arm_wall_thickness,arm_length,arm_taper_ratio,wall_length,axle_length,material_type,top_opt_thickness,L1,clevis_edge_thickness,veloMagMax,veloAngle,eigenVal1,maxMises,mass\n")

for keyVal in ['tag']: #Use taguchi only
    keyOutputDataList = []
    keyList = dataDict[keyVal]
    jobNumber = 1
    for entry in keyList[:]:
        arm_base_width,arm_base_height,arm_wall_thickness,arm_length,arm_taper_ratio,wall_length,axle_length,material_type,top_opt_thickness,L1,clevis_edge_thickness = entry
        
        material_type = int(round(material_type)) #Bucket material value index into discrete values
        
        isGood = False
        try:
            veloMagMax,veloAngle,eigenVal1,maxMises,mass = evalModel(arm_base_width,arm_base_height,arm_wall_thickness,arm_length,arm_taper_ratio,wall_length,axle_length,material_type,top_opt_thickness,L1,clevis_edge_thickness,jobNumber)
            print("Entry List '"+str(entry)+"' SUCCEEDED.")
            isGood = True
        except Exception as e:
            print("EXCEPTION: "+str(e))
            veloMagMax=-1e20
            veloAngle=0.0
            eigenVal1=-1e20
            maxMises=1e20
            mass=1e20
            print("Entry List '"+str(entry)+"' FAILED.")
        indivList = [arm_base_width,arm_base_height,arm_wall_thickness,arm_length,arm_taper_ratio,wall_length,axle_length,material_type,top_opt_thickness,L1,clevis_edge_thickness,veloMagMax,veloAngle,eigenVal1,maxMises,mass]
        keyOutputDataList.append(indivList)
    
        with open("goodBad.txt", "a") as outFile:
            outFile.write(("GOOD" if isGood else "BAD" )+"\n")
            
        with open("safetyNet.txt", "a") as outFile:
            outFile.write(",".join([str(x) for x in indivList])+"\n")
            
        jobNumber += 1
    outputDict[keyVal] = keyOutputDataList
    
#output data to JSON
with open(outputFileName, "w") as outFile:
    json.dump(outputDict,outFile)
    
print("Successfully wrote output to "+str(outputFileName))
    

