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

for keyVal in ['tag']: #Use taguchi only
    keyOutputDataList = []
    keyList = dataDict[keyVal]
    for entry in keyList[:1]:
        arm_base_width,arm_base_height,arm_wall_thickness,arm_length,arm_taper_ratio,wall_length,axle_length,material_type,top_opt_thickness,L1,L2,clevis_edge_thickness,xbeam_edge_thickness = entry
        
        material_type = int(round(material_type)) #Bucket material value index into discrete values
        
        mass,maxMises,velocity_max, angle = evalModel(arm_base_width,arm_base_height,arm_wall_thickness,arm_length,arm_taper_ratio,wall_length,axle_length,material_type,top_opt_thickness,L1,L2,clevis_edge_thickness,xbeam_edge_thickness)
        
        indivList = [arm_base_width,arm_base_height,arm_wall_thickness,arm_length,arm_taper_ratio,wall_length,axle_length,material_type,top_opt_thickness,L1,L2,clevis_edge_thickness,xbeam_edge_thickness,mass,maxMises,velocity_max,angle]
        keyOutputDataList.append(indivList)
    
    outputDict[keyVal] = keyOutputDataList
    
#output data to JSON
with open(outputFileName, "w") as outFile:
    json.dump(outputDict,outFile)
    
print("Successfully wrote output to "+str(outputFileName))
    

