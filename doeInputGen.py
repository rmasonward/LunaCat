import json
import DOEmethods as doe

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


boundsList = [[arm_base_width_min, arm_base_width_max], [arm_base_height_min, arm_base_height_max], [arm_wall_thickness_min, arm_wall_thickness_max], [arm_length_min, arm_length_max], 
[arm_taper_ratio_min, arm_taper_ratio_max], [wall_length_min, wall_length_max], [axle_length_min, axle_length_max], [material_type_min, material_type_max], [top_opt_thickness_min, top_opt_thickness_max], 
[L1_percent_min, L1_percent_max],[clevis_edge_thickness_min, clevis_edge_thickness_max]]


def createTaguchiInputFile():
    unMappedParamList = doe.Taguchi(11,"L12_Taguchi.txt")
    
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
