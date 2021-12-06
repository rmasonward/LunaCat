import json
import DOEmethods as doe

arm_base_width_min = 0.05 #m 
arm_base_width_max = 0.2 #m

arm_base_height_min = 0.01 #m
arm_base_height_max = 0.1 #m

arm_wall_thickness_min = 0.001 #m
arm_wall_thickness_max = 0.01 #m

arm_length_min = 2.5 #m
arm_length_max = 4 #m

arm_taper_ratio_min = 1
arm_taper_ratio_max = 4

wall_length_min = 1 #m
wall_length_max = 4 #m

axle_length_min = 0.3 #m
axle_length_max = 1.5 #m

#KDVs
top_opt_thickness_min = 0.05 #m
top_opt_thickness_max = 0.1 #m

L1_percent_min = 5
L1_percent_max = 35

L2_percent_min = 5
L2_percent_max = 20

clevis_edge_thickness_min = 0.01
clevis_edge_thickness_max = 0.1

xbeam_edge_thickness_min = 0.01
xbeam_edge_thickness_max = 0.1

#Axle/arm material (discrete DV)
material_type_min = 1
material_type_max = 3 


boundsList = [[arm_base_width_min, arm_base_width_max], [arm_base_height_min, arm_base_height_max], [arm_wall_thickness_min, arm_wall_thickness_max], [arm_length_min, arm_length_max], 
[arm_taper_ratio_min, arm_taper_ratio_max], [wall_length_min, wall_length_max], [axle_length_min, axle_length_max], [material_type_min, material_type_max], [top_opt_thickness_min, top_opt_thickness_max], 
[L1_percent_min, L1_percent_max], [L2_percent_min, L2_percent_max],[clevis_edge_thickness_min, clevis_edge_thickness_max], [xbeam_edge_thickness_min, xbeam_edge_thickness_max]]


def createTaguchiInputFile():
    unMappedParamList = doe.Taguchi(13,"L27_Taguchi.txt")
    
    testParamList = []
    for normedList in unMappedParamList:
        listToAdd = []
        for i,boundTuple in enumerate(boundsList):
            percentVal = normedList[i]
            mappedVal = (boundTuple[1]-boundTuple[0])*percentVal + boundTuple[0]
            listToAdd.append(mappedVal)
        testParamList.append(listToAdd)
    
    return testParamList

tagTestData = createTaguchiInputFile()

fileName = "input.json"
with open("CONFIG.txt", "r") as configFile:
    for line in configFile:
        if "INPUT FILENAME" in line:
            fileName = line.split(":")[-1].strip()

with open(fileName, "w") as outFile:
    dataToOutput = {"tag":tagTestData}
    json.dump(dataToOutput, outFile)

print("Done.")
