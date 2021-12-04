"""
Description:
This script is used to model LunaCat's final project model of a catapult.

Script by:
Peter Myers, Mason Ward, Reece Dawson, Daniel Kirby
Dept. of Aerospace Engineering
Texas A&M University
AERO 405/604 - Hartl
November 11th, 2021
"""


from abaqus import *
from abaqusConstants import *
import __main__
import section
import odbSection
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import job
import sketch
import visualization
import xyPlot
import connectorBehavior
import displayGroupOdbToolset as dgo
from math import atan2, atan, sin, cos, tan, sqrt
#from Post_P_Script import getResults



######################################
# Variable and Fixed Design Parameters
######################################

# MAKE SURE ALL VARIABLES ARE FLOATS
YOUNG_MOD = 70e9 #N/m^2
POISSON = 0.3
DENS = 2720. #kg/m^3

BASE_WIDTH = 0.14 #m ###DESIGN VARIABLE
BASE_HEIGHT = 0.06 #m ###DESIGN VARIABLE

TIP_HEIGHT = 0.02 #m

THICKNESS = 0.005 #m ###DESIGN VARIABLE
if THICKNESS*2 > min(BASE_HEIGHT,BASE_WIDTH,TIP_HEIGHT):
    raise ValueError("Thickness is too big")

BASE_LENGTH = 0.2 #m
ARM_LENGTH = 2.4 #m ###DESIGN VARIABLE
SPOON_LENGTH = 0.2 #m
SPOON_FLANGE_THICKNESS = 0.002 #m
SPOON_FLANGE_HEIGHT = 0.01 #m

WALL_THICKNESS = 0.1#m
WALL_LENGTH = 3.0 #m
WALL_HEIGHT = 0.4 #m


# Axle Beam Variables
S1 = 0.1 #m  used to sketch the cross section
AXLE_LENGTH = .6 #m ###DESIGN VARIABLE
#This breaks when changing value over 1
AXLE_DIAMETER = 0.2 #m  this is defined as the hypotenuse of the square cross section


# Cross Member Variables
WALL_SEP_LENGTH = AXLE_LENGTH-2*WALL_THICKNESS #m

XBEAM_HEIGHT = 0.7 #m
XBEAM_FLAT_LENGTH = 0.02 #m
XBEAM_THICKNESS = 0.1 #m
XBEAM_AXLE_DIA = 0.2 #m
OFFSET = sqrt((BASE_HEIGHT**2)/8) ###DESIGN VARIABLE
CLEVIS_RAD = XBEAM_AXLE_DIA/2
Clev_Opp = CLEVIS_RAD * cos(pi/4)


# New variables from top opt
#### 5 design variables for top opt
ToptThickness = 0.05
L1 = 0.5
L2 = 0.3
CLEVIS_EDGE_THICK = 0.05
XBEAM_EDGE_THICK = 0.05

#fixed variables
CLEVIS_DIST = (-WALL_LENGTH/2)+CLEVIS_RAD+CLEVIS_EDGE_THICK   ###MIDDLE OF CIRCULAR NOTCH
WALL_SIDE_X = (WALL_LENGTH/2)- XBEAM_EDGE_THICK - (AXLE_DIAMETER/2) ###MIDDLE OF SQUARE NOTCH
###variables for square notch
P1_X = WALL_SIDE_X #m
P1_Y = -AXLE_DIAMETER/2 #m
P2_X = WALL_SIDE_X + AXLE_DIAMETER/2
P2_Y = 0 #m
P3_X = WALL_SIDE_X #m
P3_Y = AXLE_DIAMETER/2 #m
P4_X = WALL_SIDE_X - AXLE_DIAMETER/2 #m
P4_Y = 0 #m
## Variables for top opt cut
X1 = WALL_LENGTH/2-AXLE_DIAMETER- 2*XBEAM_EDGE_THICK
Y1 = -WALL_HEIGHT/2 + ToptThickness
X2 = X1 - L1
Y2 = Y1
X3 = X2
Y3 = 0
X4 = X3 - L2
Y4 = WALL_HEIGHT/2 - ToptThickness
X5 = X1
Y5 = Y4
X6 = X1
Y6 = 0


#Part Names
OUTER_ARM_NAME = "OUTER_ARM_SOLID"
INNER_ARM_NAME = "ARM_VOID_SPACE"
ARM_NAME = "ARM"
XBEAM_NAME = "CROSSMEMBER"
##########################
# FEA Modeling Parameters
# (e.g., mesh seeds, step times, etc)
##########################

####################################
### Calculated Properties/Values ###
####################################

session.journalOptions.setValues(replayGeometry=COORDINATE,recoverGeometry=COORDINATE)
session.viewports['Viewport: 1'].assemblyDisplay.geometryOptions.setValues(
    datumPlanes=OFF)

### Write data file column headings
DataFile = open('PostData.txt','w')
DataFile.write('seedSize tipDisp \n')
DataFile.close()

###############################
### Generation of FEA Model ###
###############################

### Note: If you create a loop, START it here 
### (i.e., the contents of the loop below are intended)


### Scripting the entire model allows its entire
### contents to be packaged into this single file.
Mdb() 


##################################  
# Sketch Geometry and Create Parts
##################################
print('Sketching/Creating the part')

###########################
### TopOpted SideWall 1 ###
###########################
s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
    sheetSize=10.0)
g, v, d2, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
s1.setPrimaryObject(option=STANDALONE)

s1.Line(point1=(-WALL_LENGTH/5, WALL_HEIGHT/2), point2=(CLEVIS_DIST, CLEVIS_RAD + CLEVIS_EDGE_THICK))
s1.Line(point1=(-WALL_LENGTH/5, -WALL_HEIGHT/2), point2=(CLEVIS_DIST, -(CLEVIS_RAD + CLEVIS_EDGE_THICK)))
s1.Arc3Points(point1=(CLEVIS_DIST, CLEVIS_RAD+CLEVIS_EDGE_THICK), point2=(CLEVIS_DIST, -(CLEVIS_RAD+CLEVIS_EDGE_THICK)),
    point3=(-WALL_LENGTH/2, 0.0))
s1.Line(point1=(-WALL_LENGTH/5, -WALL_HEIGHT/2), point2=(WALL_SIDE_X, -(WALL_HEIGHT/2)))
s1.Line(point1=(WALL_SIDE_X, -(WALL_HEIGHT/2)), point2=(WALL_LENGTH/2, 0.0))
s1.Line(point1=(WALL_LENGTH/2, 0.0), point2=(WALL_SIDE_X, WALL_HEIGHT/2))
s1.Line(point1=(WALL_SIDE_X, WALL_HEIGHT/2), point2=(-WALL_LENGTH/5, WALL_HEIGHT/2))
p = mdb.models['Model-1'].Part(name='Side_wall_1', dimensionality=THREE_D, 
    type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Side_wall_1']
p.BaseSolidExtrude(sketch=s1, depth=WALL_THICKNESS)
s1.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['Side_wall_1']

#Datum plane on side wall 1
p = mdb.models['Model-1'].parts['Side_wall_1']
p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=0.0)

#Partitioning 
p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=WALL_SIDE_X)
c = p.cells
pickedCells = c.findAt(((0.0, 0.0, 0.0), ))
d1 = p.datums
p.PartitionCellByDatumPlane(datumPlane=d1[3], cells=pickedCells)
p = mdb.models['Model-1'].parts['Side_wall_1']
p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=WALL_THICKNESS/2)
c = p.cells
pickedCells = c.findAt(((0.0, 0.0, 0.0), ), ((WALL_LENGTH/2-AXLE_DIAMETER/4, 0.0, 0.0), ))
d2 = p.datums
p.PartitionCellByDatumPlane(datumPlane=d2[5], cells=pickedCells)
p = mdb.models['Model-1'].parts['Side_wall_1']
mdb.models['Model-1'].sketches.changeKey(fromName='__profile__', 
    toName='__save__')
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=10.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
s.retrieveSketch(sketch=mdb.models['Model-1'].sketches['__save__'])
s.Line(point1=(P1_X, P1_Y), point2=(P2_X, P2_Y))
s.Line(point1=(P2_X, P2_Y), point2=(P3_X, P3_Y))
s.Line(point1=(P3_X, P3_Y), point2=(P4_X, P4_Y))
s.Line(point1=(P4_X, P4_Y), point2=(P1_X, P1_Y))
p = mdb.models['Model-1'].Part(name='SquareCutOut', dimensionality=THREE_D, 
    type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['SquareCutOut']
p.BaseSolidExtrude(sketch=s, depth=WALL_THICKNESS)
s.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['SquareCutOut']

## Clevis Cut in Side Wall 1
mdb.models['Model-1'].sketches.changeKey(fromName='__profile__', 
    toName='__save__')
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=10.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
s.retrieveSketch(sketch=mdb.models['Model-1'].sketches['__save__'])
del mdb.models['Model-1'].sketches['__save__']
s.CircleByCenterPerimeter(center=(CLEVIS_DIST, 0.0), point1=(CLEVIS_DIST+CLEVIS_RAD, 
    0.0))
p = mdb.models['Model-1'].Part(name='XbeamCutOut', dimensionality=THREE_D, 
    type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['XbeamCutOut']
p.BaseSolidExtrude(sketch=s, depth=WALL_THICKNESS)
s.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['XbeamCutOut']
#Creating top opt cut 
mdb.models['Model-1'].sketches.changeKey(fromName='__profile__', 
    toName='__save__')
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=10.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
s.retrieveSketch(sketch=mdb.models['Model-1'].sketches['__save__'])
del mdb.models['Model-1'].sketches['__save__']
s.Line(point1=(X1, Y1), point2=(X2, Y2))
s.Line(point1=(X2, Y2), point2=(X3,Y3))
s.Line(point1=(X3,Y3), point2=(X4,Y4))
s.Line(point1=(X4,Y4), point2=(X5,Y5))
s.Line(point1=(X5,Y5), point2=(X6,Y6))
s.Line(point1=(X6,Y6), point2=(X1, Y1))
p = mdb.models['Model-1'].Part(name='TopOptCut', dimensionality=THREE_D, 
    type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['TopOptCut']
p.BaseSolidExtrude(sketch=s, depth=WALL_THICKNESS)
s.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['TopOptCut']

## Cleaning up model
p1 = mdb.models['Model-1'].parts['XbeamCutOut']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
del mdb.models['Model-1'].parts['XbeamCutOut']
p = mdb.models['Model-1'].parts['Side_wall_1']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
p1 = mdb.models['Model-1'].parts['SquareCutOut']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
del mdb.models['Model-1'].parts['SquareCutOut']
p = mdb.models['Model-1'].parts['Side_wall_1']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
del mdb.models['Model-1'].parts['Side_wall_1']
p = mdb.models['Model-1'].parts['TopOptCut']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
mdb.models['Model-1'].parts.changeKey(fromName='TopOptCut', 
    toName='Side_wall_1')

########################################
##Arm
########################################
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
    sheetSize=200.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
s.Line(point1=(0.0, BASE_HEIGHT/2), point2=(0.0, -BASE_HEIGHT/2))
s.Line(point1=(0.0, -BASE_HEIGHT/2), point2=(BASE_LENGTH, -BASE_HEIGHT/2))
s.Line(point1=(BASE_LENGTH, -BASE_HEIGHT/2), point2=(BASE_LENGTH+ARM_LENGTH, -TIP_HEIGHT/2))
s.Line(point1=(BASE_LENGTH+ARM_LENGTH, -TIP_HEIGHT/2), point2=(BASE_LENGTH+ARM_LENGTH, TIP_HEIGHT/2))
s.Line(point1=(BASE_LENGTH+ARM_LENGTH, TIP_HEIGHT/2), point2=(BASE_LENGTH, BASE_HEIGHT/2))
s.Line(point1=(BASE_LENGTH, BASE_HEIGHT/2), point2=(0.0, BASE_HEIGHT/2))


p = mdb.models['Model-1'].Part(name=OUTER_ARM_NAME, dimensionality=THREE_D, 
    type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts[OUTER_ARM_NAME]
p.BaseSolidExtrude(sketch=s, depth=BASE_WIDTH)
s.unsetPrimaryObject()
del mdb.models['Model-1'].sketches['__profile__']

#Now make the inner arm space that will be subtracted (the inner void volume)
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
    sheetSize=200.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
s.Line(point1=(0.0, BASE_HEIGHT/2-THICKNESS), point2=(0.0, -(BASE_HEIGHT/2-THICKNESS)))
s.Line(point1=(0.0, -(BASE_HEIGHT/2-THICKNESS)), point2=(BASE_LENGTH, -(BASE_HEIGHT/2-THICKNESS)))
s.Line(point1=(BASE_LENGTH, -(BASE_HEIGHT/2-THICKNESS)), point2=(BASE_LENGTH+ARM_LENGTH, -(TIP_HEIGHT/2-THICKNESS)))
s.Line(point1=(BASE_LENGTH+ARM_LENGTH, -(TIP_HEIGHT/2-THICKNESS)), point2=(BASE_LENGTH+ARM_LENGTH, TIP_HEIGHT/2-THICKNESS))
s.Line(point1=(BASE_LENGTH+ARM_LENGTH, TIP_HEIGHT/2-THICKNESS), point2=(BASE_LENGTH, BASE_HEIGHT/2-THICKNESS))
s.Line(point1=(BASE_LENGTH, BASE_HEIGHT/2-THICKNESS), point2=(0.0, BASE_HEIGHT/2-THICKNESS))


p = mdb.models['Model-1'].Part(name=INNER_ARM_NAME, dimensionality=THREE_D, 
    type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts[INNER_ARM_NAME]
p.BaseSolidExtrude(sketch=s, depth=(BASE_WIDTH-2*THICKNESS))
s.unsetPrimaryObject()
del mdb.models['Model-1'].sketches['__profile__']

#Use boolean cut to make the final arm here
OUTER_ARM_INST_NAME = OUTER_ARM_NAME+'-1'
INNER_ARM_INST_NAME = INNER_ARM_NAME+'-1'

a = mdb.models['Model-1'].rootAssembly
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
p = mdb.models['Model-1'].parts[INNER_ARM_NAME]
a.Instance(name=INNER_ARM_INST_NAME, part=p, dependent=ON)

p = mdb.models['Model-1'].parts[OUTER_ARM_NAME]
a.Instance(name=OUTER_ARM_INST_NAME, part=p, dependent=ON)

a.translate(instanceList=(INNER_ARM_INST_NAME, ), vector=(0.0, 0.0, THICKNESS))

#Cut out the inner void space to make the final arm
a.InstanceFromBooleanCut(name=ARM_NAME, 
    instanceToBeCut=mdb.models['Model-1'].rootAssembly.instances[OUTER_ARM_INST_NAME], 
    cuttingInstances=(a.instances[INNER_ARM_INST_NAME], ), 
    originalInstances=SUPPRESS)

#make the "spoon" section
p = mdb.models['Model-1'].parts[ARM_NAME]
f, e = p.faces, p.edges
t = p.MakeSketchTransform(sketchPlane=f.findAt(coordinates=(BASE_LENGTH+ARM_LENGTH, 0.0, 
    THICKNESS/2)), sketchUpEdge=e.findAt(coordinates=(BASE_LENGTH+ARM_LENGTH, 0.0, 0.0)), 
    sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(BASE_LENGTH+ARM_LENGTH, 0.0, BASE_WIDTH/2))
s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
    sheetSize=0.282, gridSpacing=0.007, transform=t)
g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
s1.sketchOptions.setValues(decimalPlaces=3)
s1.setPrimaryObject(option=SUPERIMPOSE)
p = mdb.models['Model-1'].parts[ARM_NAME]
p.projectReferencesOntoSketch(sketch=s1, filter=COPLANAR_EDGES)
s1.rectangle(point1=(-BASE_WIDTH/2, TIP_HEIGHT/2), point2=(BASE_WIDTH/2, -TIP_HEIGHT/2))
f1, e1 = p.faces, p.edges
p.SolidExtrude(sketchPlane=f1.findAt(coordinates=(BASE_LENGTH+ARM_LENGTH, 0.0, THICKNESS/2)), 
    sketchUpEdge=e1.findAt(coordinates=(BASE_LENGTH+ARM_LENGTH, 0.0, 0.0)), 
    sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, sketch=s1, depth=SPOON_LENGTH, 
    flipExtrudeDirection=OFF)
s1.unsetPrimaryObject()
del mdb.models['Model-1'].sketches['__profile__']

#Make the payload holders here
t = p.MakeSketchTransform(sketchPlane=f.findAt(coordinates=(ARM_LENGTH+BASE_LENGTH-SPOON_LENGTH/2, 0.0, 
    BASE_WIDTH)), sketchUpEdge=e.findAt(coordinates=(ARM_LENGTH+BASE_LENGTH+SPOON_LENGTH, 0.0, BASE_WIDTH)), 
    sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(ARM_LENGTH+BASE_LENGTH, TIP_HEIGHT/2, 
    BASE_WIDTH))
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=11.2, 
    gridSpacing=0.28, transform=t)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.sketchOptions.setValues(decimalPlaces=3)
s.setPrimaryObject(option=SUPERIMPOSE)
p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)

#proximal flange
s.Line(point1=(0.0,0.0), point2=(0.0, SPOON_FLANGE_HEIGHT))
s.Line(point1=(0.0, SPOON_FLANGE_HEIGHT), point2=(SPOON_FLANGE_THICKNESS, SPOON_FLANGE_HEIGHT))
s.Line(point1=(SPOON_FLANGE_THICKNESS, SPOON_FLANGE_HEIGHT), point2=(SPOON_FLANGE_THICKNESS, SPOON_FLANGE_THICKNESS))
s.Line(point1=(SPOON_FLANGE_THICKNESS, SPOON_FLANGE_THICKNESS), point2=(SPOON_FLANGE_HEIGHT, SPOON_FLANGE_THICKNESS))
s.Line(point1=(SPOON_FLANGE_HEIGHT, SPOON_FLANGE_THICKNESS), point2=(SPOON_FLANGE_HEIGHT, 0.0))
s.Line(point1=(SPOON_FLANGE_HEIGHT, 0.0), point2=(0.0, 0.0))

#distal flange
s.Line(point1=(SPOON_LENGTH-SPOON_FLANGE_HEIGHT,0.0), point2=(SPOON_LENGTH-SPOON_FLANGE_HEIGHT,SPOON_FLANGE_THICKNESS))
s.Line(point1=(SPOON_LENGTH-SPOON_FLANGE_HEIGHT,SPOON_FLANGE_THICKNESS), point2=(SPOON_LENGTH-SPOON_FLANGE_THICKNESS,SPOON_FLANGE_THICKNESS))
s.Line(point1=(SPOON_LENGTH-SPOON_FLANGE_THICKNESS,SPOON_FLANGE_THICKNESS), point2=(SPOON_LENGTH-SPOON_FLANGE_THICKNESS,SPOON_FLANGE_HEIGHT))
s.Line(point1=(SPOON_LENGTH-SPOON_FLANGE_THICKNESS,SPOON_FLANGE_HEIGHT), point2=(SPOON_LENGTH,SPOON_FLANGE_HEIGHT))
s.Line(point1=(SPOON_LENGTH,SPOON_FLANGE_HEIGHT), point2=(SPOON_LENGTH,0.0))
s.Line(point1=(SPOON_LENGTH,0.0), point2=(SPOON_LENGTH-SPOON_FLANGE_HEIGHT,0.0))

p.SolidExtrude(sketchPlane=f1.findAt(coordinates=(ARM_LENGTH+BASE_LENGTH-SPOON_LENGTH/2, 0.0, BASE_WIDTH)), 
    sketchUpEdge=e1.findAt(coordinates=(ARM_LENGTH+BASE_LENGTH+SPOON_LENGTH, 0.0, BASE_WIDTH)), 
    sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, sketch=s, depth=BASE_WIDTH, 
    flipExtrudeDirection=ON)
s.unsetPrimaryObject()
del mdb.models['Model-1'].sketches['__profile__']
########################################################################

##Connecting beam
s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
    sheetSize=10.0)
g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
s1.setPrimaryObject(option=STANDALONE)
s1.Line(point1=(0.0, 0.0), point2=(AXLE_DIAMETER/2, AXLE_DIAMETER/2))
s1.Line(point1=(AXLE_DIAMETER/2, AXLE_DIAMETER/2), point2=(0.0, AXLE_DIAMETER))
# s1.PerpendicularConstraint(entity1=g.findAt((0.1, 0.1)), entity2=g.findAt((0.1, 
    # 0.3)), addUndoState=False)
s1.Line(point1=(0.0, AXLE_DIAMETER), point2=(-AXLE_DIAMETER/2, AXLE_DIAMETER/2))
# s1.PerpendicularConstraint(entity1=g.findAt((0.1, 0.3)), entity2=g.findAt((
    # -0.1, 0.3)), addUndoState=False)
s1.Line(point1=(-AXLE_DIAMETER/2, AXLE_DIAMETER/2), point2=(0.0, 0.0))
# s1.PerpendicularConstraint(entity1=g.findAt((-0.1, 0.3)), entity2=g.findAt((
    # -0.1, 0.1)), addUndoState=False)
p = mdb.models['Model-1'].Part(name='Connector_beam', dimensionality=THREE_D, 
    type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Connector_beam']
p.BaseSolidExtrude(sketch=s1, depth=AXLE_LENGTH)
s1.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['Connector_beam']
del mdb.models['Model-1'].sketches['__profile__']
p = mdb.models['Model-1'].parts['Connector_beam']
v2 = p.vertices

#Datum plane
p.DatumPlaneByThreePoints(point1=v2.findAt(coordinates=(0.0, AXLE_DIAMETER, AXLE_LENGTH)), 
    point2=v2.findAt(coordinates=(AXLE_DIAMETER/2, AXLE_DIAMETER/2, AXLE_LENGTH)), point3=v2.findAt(
    coordinates=(AXLE_DIAMETER/2, AXLE_DIAMETER/2, 0.0)))

#Extruded slot in connector beam
p = mdb.models['Model-1'].parts['Connector_beam']
e, d = p.edges, p.datums
t = p.MakeSketchTransform(sketchPlane=d[2], sketchUpEdge=e.findAt(coordinates=(
    AXLE_DIAMETER/2, AXLE_DIAMETER/2, AXLE_LENGTH/2)), sketchPlaneSide=SIDE1, sketchOrientation=BOTTOM, origin=(
    AXLE_DIAMETER/4, 1.5*AXLE_DIAMETER/2, AXLE_LENGTH/2))
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=4.46, 
    gridSpacing=0.11, transform=t)
g, v, d1, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=SUPERIMPOSE)
p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
s.rectangle(point1=(-BASE_WIDTH/2, BASE_HEIGHT/2), point2=(BASE_WIDTH/2, -BASE_HEIGHT/2))
e1, d2 = p.edges, p.datums
p.CutExtrude(sketchPlane=d2[2], sketchUpEdge=e1.findAt(coordinates=(AXLE_DIAMETER/2, AXLE_DIAMETER/2, 
    AXLE_LENGTH/2)), sketchPlaneSide=SIDE1, sketchOrientation=BOTTOM, sketch=s, 
    flipExtrudeDirection=OFF)
s.unsetPrimaryObject()
del mdb.models['Model-1'].sketches['__profile__']


    
##Clevis Triangle
 #define vertical sketch offset
vertOffset = XBEAM_AXLE_DIA*0.75
slope = (XBEAM_HEIGHT-vertOffset)/(WALL_SEP_LENGTH/2-XBEAM_FLAT_LENGTH/2)
insetDist = WALL_SEP_LENGTH/2 - (XBEAM_HEIGHT-1.75*vertOffset)/slope

s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
    sheetSize=200.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)

print("Axle length: "+str(AXLE_LENGTH))
print("Wall sep dist: "+str(WALL_SEP_LENGTH))
print((AXLE_LENGTH-WALL_SEP_LENGTH)/2)

s.Line(point1=((AXLE_LENGTH-WALL_SEP_LENGTH)/2, 0.0), point2=((AXLE_LENGTH-WALL_SEP_LENGTH)/2+WALL_SEP_LENGTH, 0.0))
s.Line(point1=((AXLE_LENGTH-WALL_SEP_LENGTH)/2+WALL_SEP_LENGTH, 0.0), point2=((AXLE_LENGTH-WALL_SEP_LENGTH)/2+WALL_SEP_LENGTH, vertOffset))
s.Line(point1=((AXLE_LENGTH-WALL_SEP_LENGTH)/2+WALL_SEP_LENGTH, vertOffset), point2=((AXLE_LENGTH-WALL_SEP_LENGTH)/2+WALL_SEP_LENGTH/2+XBEAM_FLAT_LENGTH/2, XBEAM_HEIGHT))
s.Line(point1=((AXLE_LENGTH-WALL_SEP_LENGTH)/2+WALL_SEP_LENGTH/2+XBEAM_FLAT_LENGTH/2, XBEAM_HEIGHT), point2=((AXLE_LENGTH-WALL_SEP_LENGTH)/2+WALL_SEP_LENGTH/2-XBEAM_FLAT_LENGTH/2, XBEAM_HEIGHT))
s.Line(point1=((AXLE_LENGTH-WALL_SEP_LENGTH)/2+WALL_SEP_LENGTH/2-XBEAM_FLAT_LENGTH/2, XBEAM_HEIGHT), point2=((AXLE_LENGTH-WALL_SEP_LENGTH)/2, vertOffset))
s.Line(point1=((AXLE_LENGTH-WALL_SEP_LENGTH)/2, vertOffset), point2=((AXLE_LENGTH-WALL_SEP_LENGTH)/2, 0.0))

#make cutout
s.Line(point1=((AXLE_LENGTH-WALL_SEP_LENGTH)/2+WALL_SEP_LENGTH/2, XBEAM_HEIGHT-0.75*vertOffset), point2=((AXLE_LENGTH-WALL_SEP_LENGTH)/2+insetDist, vertOffset))
s.Line(point1=((AXLE_LENGTH-WALL_SEP_LENGTH)/2+insetDist, vertOffset), point2=((AXLE_LENGTH-WALL_SEP_LENGTH)/2+WALL_SEP_LENGTH-insetDist, vertOffset))
s.Line(point1=((AXLE_LENGTH-WALL_SEP_LENGTH)/2+WALL_SEP_LENGTH-insetDist, vertOffset), point2=((AXLE_LENGTH-WALL_SEP_LENGTH)/2+WALL_SEP_LENGTH/2, XBEAM_HEIGHT-0.75*vertOffset))

p = mdb.models['Model-1'].Part(name=XBEAM_NAME, dimensionality=THREE_D, 
    type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts[XBEAM_NAME]
p.BaseSolidExtrude(sketch=s, depth=XBEAM_THICKNESS)
s.unsetPrimaryObject()
del mdb.models['Model-1'].sketches['__profile__']

#extrude axle
p = mdb.models['Model-1'].parts[XBEAM_NAME]
f, e = p.faces, p.edges
t = p.MakeSketchTransform(sketchPlane=f.findAt(coordinates=((AXLE_LENGTH-WALL_SEP_LENGTH)/2, vertOffset/2, 
    XBEAM_THICKNESS/2)), sketchUpEdge=e.findAt(coordinates=((AXLE_LENGTH-WALL_SEP_LENGTH)/2, vertOffset/2, XBEAM_THICKNESS)), 
    sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=((AXLE_LENGTH-WALL_SEP_LENGTH)/2, vertOffset/2, XBEAM_THICKNESS/2))
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=1.31, 
    gridSpacing=0.03, transform=t)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=SUPERIMPOSE)
p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)
s.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(XBEAM_AXLE_DIA/2, 0.0))
f1, e1 = p.faces, p.edges
p.SolidExtrude(sketchPlane=f1.findAt(coordinates=((AXLE_LENGTH-WALL_SEP_LENGTH)/2, vertOffset/2, XBEAM_THICKNESS/2)), 
    sketchUpEdge=e1.findAt(coordinates=((AXLE_LENGTH-WALL_SEP_LENGTH)/2, vertOffset/2, XBEAM_THICKNESS)), 
    sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, sketch=s, depth=WALL_SEP_LENGTH+(AXLE_LENGTH-WALL_SEP_LENGTH)/2, 
    flipExtrudeDirection=ON)
s.unsetPrimaryObject()
del mdb.models['Model-1'].sketches['__profile__']

#extrude last bit of axle
t = p.MakeSketchTransform(sketchPlane=f.findAt(coordinates=((AXLE_LENGTH-WALL_SEP_LENGTH)/2, vertOffset/2, 
    XBEAM_THICKNESS/2)), sketchUpEdge=e.findAt(coordinates=((AXLE_LENGTH-WALL_SEP_LENGTH)/2+WALL_SEP_LENGTH/2-XBEAM_FLAT_LENGTH/2, XBEAM_HEIGHT, XBEAM_THICKNESS/2)), 
    sketchPlaneSide=SIDE1, sketchOrientation=TOP, origin=((AXLE_LENGTH-WALL_SEP_LENGTH)/2, vertOffset/2, XBEAM_THICKNESS/2))
s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
    sheetSize=1.54, gridSpacing=0.03, transform=t)
g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
s1.setPrimaryObject(option=SUPERIMPOSE)
p.projectReferencesOntoSketch(sketch=s1, filter=COPLANAR_EDGES)
s1.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(XBEAM_AXLE_DIA/2, 0.0))
f1, e1 = p.faces, p.edges
p.SolidExtrude(sketchPlane=f1.findAt(coordinates=((AXLE_LENGTH-WALL_SEP_LENGTH)/2, vertOffset/2, XBEAM_THICKNESS/2)),
    sketchUpEdge=e1.findAt(coordinates=((AXLE_LENGTH-WALL_SEP_LENGTH)/2+WALL_SEP_LENGTH/2-XBEAM_FLAT_LENGTH/2, XBEAM_HEIGHT, XBEAM_THICKNESS/2)), 
    sketchPlaneSide=SIDE1, sketchOrientation=TOP, sketch=s1, depth=(AXLE_LENGTH-WALL_SEP_LENGTH)/2, 
    flipExtrudeDirection=OFF)
s1.unsetPrimaryObject()
del mdb.models['Model-1'].sketches['__profile__']

#add fillets
# e = p.edges
p.Round(radius=AXLE_LENGTH/8, edgeList=(e.findAt(coordinates=((AXLE_LENGTH-WALL_SEP_LENGTH)/2+WALL_SEP_LENGTH/2, XBEAM_HEIGHT-vertOffset*0.75, XBEAM_THICKNESS/2)), 
     ))   
# #create reference point
# p.DatumPointByCoordinate(coords=((AXLE_LENGTH-WALL_SEP_LENGTH)/2+WALL_SEP_LENGTH/2, XBEAM_HEIGHT, XBEAM_THICKNESS/2))

# #now build partitions
# c = p.cells
# pickedCells = c.findAt((((AXLE_LENGTH-WALL_SEP_LENGTH)/2+WALL_SEP_LENGTH/2, XBEAM_HEIGHT, XBEAM_THICKNESS/2), ))
# e, v, d = p.edges, p.vertices, p.datums
# p.PartitionCellByPlanePointNormal(point=d[5], normal=e.findAt(coordinates=((AXLE_LENGTH-WALL_SEP_LENGTH)/2+WALL_SEP_LENGTH/2-XBEAM_FLAT_LENGTH/2, XBEAM_HEIGHT, XBEAM_THICKNESS/4)), cells=pickedCells)
# # #next partition
# c = p.cells
# pickedCells = c.findAt((((AXLE_LENGTH-WALL_SEP_LENGTH)/2+WALL_SEP_LENGTH/2, XBEAM_HEIGHT, XBEAM_THICKNESS/4), ), (((AXLE_LENGTH-WALL_SEP_LENGTH)/2+WALL_SEP_LENGTH/2, XBEAM_HEIGHT, 3*XBEAM_THICKNESS/4), ))
# p.PartitionCellByPlanePointNormal(point=d[5], normal=e1.findAt(coordinates=((AXLE_LENGTH-WALL_SEP_LENGTH)/2+WALL_SEP_LENGTH/2+XBEAM_FLAT_LENGTH/4, XBEAM_HEIGHT, 0.0)), cells=pickedCells)

# #create pull point set
# v = p.vertices
# verts = v.findAt((((AXLE_LENGTH-WALL_SEP_LENGTH)/2+WALL_SEP_LENGTH/2, XBEAM_HEIGHT, XBEAM_THICKNESS/2), ))
# p.Set(vertices=verts, name='XBEAM_PULL_POINT')



##################################  
# Create Material
##################################  

print('Creating the Materials')
mdb.models['Model-1'].Material(name='Aluminum')
mdb.models['Model-1'].materials['Aluminum'].Density(table=((2710.0, ), ))
mdb.models['Model-1'].materials['Aluminum'].Elastic(table=((10000000000.0, 
    0.3), ))

mdb.models['Model-1'].Material(name='Graphite Epoxy AS/3501')
mdb.models['Model-1'].materials['Graphite Epoxy AS/3501'].Elastic(type=LAMINA, 
    table=((137894900000.0, 8963168000.0, 2068423000.0, 
    6894745000.0, 6894745000.0, 6894745000.0), ))



##################################  
#Create/Assign Section
##################################  

print('Creating the Sections')
#Creating section 1
p = mdb.models['Model-1'].parts['Side_wall_1']
mdb.models['Model-1'].HomogeneousSolidSection(name='Section-1', 
    material='Aluminum', thickness=None)
    
#Creating Section 2
mdb.models['Model-1'].HomogeneousSolidSection(name='Section-2', 
    material='Graphite Epoxy AS/3501', thickness=None)

#Assigning Section 1 to SIDE WALL 1
p = mdb.models['Model-1'].parts['Side_wall_1']
c = p.cells
cells = c.findAt(((WALL_LENGTH/2 - XBEAM_EDGE_THICK/2, 0.0, WALL_THICKNESS/2), ))
region = p.Set(cells=cells, name='Set-1')
p = mdb.models['Model-1'].parts['Side_wall_1']
p.SectionAssignment(region=region, sectionName='Section-1', offset=0.0, 
    offsetType=MIDDLE_SURFACE, offsetField='', 
    thicknessAssignment=FROM_SECTION)
    

#Assigning Section 1 to AXLE
p = mdb.models['Model-1'].parts['Connector_beam']
c = p.cells
cells = c.findAt(((0.0, 0.0, 0.0), ))
region = p.Set(cells=cells, name='Set-1')
p = mdb.models['Model-1'].parts['Connector_beam']
p.SectionAssignment(region=region, sectionName='Section-1', offset=0.0, 
    offsetType=MIDDLE_SURFACE, offsetField='', 
    thicknessAssignment=FROM_SECTION)

#Assigning Section 1 to CLEVIS
p = mdb.models['Model-1'].parts['CROSSMEMBER']
c = p.cells
cells = c.findAt(((0.0, 0.0, 0.0), ))
region = p.Set(cells=cells, name='Set-1')
p = mdb.models['Model-1'].parts['CROSSMEMBER']
p.SectionAssignment(region=region, sectionName='Section-1', offset=0.0, 
    offsetType=MIDDLE_SURFACE, offsetField='', 
    thicknessAssignment=FROM_SECTION)

#Assigning Section 2 to ARM
p = mdb.models['Model-1'].parts['ARM']
c = p.cells
cells = c.findAt(((0.0, 0.0, 0.0), ))
region = p.Set(cells=cells, name='Set-1')
p = mdb.models['Model-1'].parts['ARM']
p.SectionAssignment(region=region, sectionName='Section-2', offset=0.0, 
    offsetType=MIDDLE_SURFACE, offsetField='', 
    thicknessAssignment=FROM_SECTION)


##################################  
# Create Composite Layup
##################################  

armDatumIndices = {} ### Creates array to reference datum plane
p = mdb.models['Model-1'].parts['ARM']
v, e = p.vertices, p.edges
d = p.datums
    
p.DatumCsysByThreePoints(point2=v.findAt(coordinates=(0.0, BASE_HEIGHT/2, 0.0)), #Datum 28
    name='Datum csys-1', coordSysType=CARTESIAN, origin=(0.0, 0.0, 0.0), 
    point1=p.InterestingPoint(edge=e.findAt(coordinates=(ARM_LENGTH+BASE_LENGTH+SPOON_LENGTH, 0.0, 0.0)), 
    rule=MIDDLE))
    
armDatumIndices['Composite-csys'] = d.keys()[-1]
layupOrientation = d[armDatumIndices['Composite-csys']]

c = p.cells
cells = c.findAt(((0.0, 0.0, 0.0), ), ((BASE_LENGTH+ARM_LENGTH+SPOON_LENGTH, 0.0, 0.0), 
    ))
region1=regionToolset.Region(cells=cells)
region2=regionToolset.Region(cells=cells)
region3=regionToolset.Region(cells=cells)
compositeLayup = mdb.models['Model-1'].parts['ARM'].CompositeLayup(
    name='CompositeLayup-1', description='', elementType=SOLID, 
    symmetric=False, thicknessAssignment=FROM_SECTION)
compositeLayup.ReferenceOrientation(orientationType=SYSTEM, 
    localCsys=layupOrientation, fieldName='', 
    additionalRotationType=ROTATION_NONE, angle=0.0, 
    additionalRotationField='', axis=AXIS_3, stackDirection=STACK_3)
compositeLayup.CompositePly(suppressed=False, plyName='Ply-1', region=region1, 
    material='Graphite Epoxy AS/3501', thicknessType=SPECIFY_THICKNESS, 
    thickness=0.1, orientationType=ANGLE_0, 
    additionalRotationType=ROTATION_NONE, additionalRotationField='', 
    axis=AXIS_3, angle=0.0, numIntPoints=3)
compositeLayup.CompositePly(suppressed=False, plyName='Ply-2', region=region2, 
    material='Graphite Epoxy AS/3501', thicknessType=SPECIFY_THICKNESS, 
    thickness=0.1, orientationType=ANGLE_0, 
    additionalRotationType=ROTATION_NONE, additionalRotationField='', 
    axis=AXIS_3, angle=45, numIntPoints=3)
compositeLayup.CompositePly(suppressed=False, plyName='Ply-3', region=region3, 
    material='Graphite Epoxy AS/3501', thicknessType=SPECIFY_THICKNESS, 
    thickness=0.1, orientationType=ANGLE_0, 
    additionalRotationType=ROTATION_NONE, additionalRotationField='', 
    axis=AXIS_3, angle=-45, numIntPoints=3)



##################################  
#Defining the face partitions
##################################  
print('Partitioning part')

#### ARM ####
#Datum Planes for Assembly
p = mdb.models['Model-1'].parts[ARM_NAME]
p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=BASE_LENGTH)
armDatumIndices['armplane1'] = p.datums.keys()[-1]
c = p.cells
pickedCells = c.findAt(((BASE_LENGTH/2, 0.0, THICKNESS/2), ))
d = p.datums
p.PartitionCellByDatumPlane(datumPlane=d[armDatumIndices['armplane1']], cells=pickedCells)
#5

#partition "spoon" section first
c = p.cells
pickedCells = c.findAt(((BASE_LENGTH+ARM_LENGTH+SPOON_LENGTH/2, 0.0, BASE_WIDTH/2), ))
p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=BASE_LENGTH+ARM_LENGTH)
armDatumIndices['armplane2'] = p.datums.keys()[-1]
p.PartitionCellByDatumPlane(datumPlane=d[armDatumIndices['armplane2']], cells=pickedCells)

#Partitioning side walls 
p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=THICKNESS)
armDatumIndices['armplane3'] = p.datums.keys()[-1]
p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=BASE_WIDTH-THICKNESS)
armDatumIndices['armplane4'] = p.datums.keys()[-1]
c = p.cells
pickedCells = c.findAt(((BASE_LENGTH/2, 0.0, THICKNESS/2), ), ((BASE_LENGTH+ARM_LENGTH/2, 0.0, THICKNESS/2), ))
p.PartitionCellByDatumPlane(datumPlane=d[armDatumIndices['armplane3']], cells=pickedCells)
c = p.cells
p = mdb.models['Model-1'].parts['ARM']
c = p.cells
pickedCells = c.findAt(((BASE_LENGTH/2, 0.0, BASE_WIDTH-THICKNESS/2), ), ((BASE_LENGTH+ARM_LENGTH/2, 0.0, BASE_WIDTH-THICKNESS/2), ))
p.PartitionCellByDatumPlane(datumPlane=d[armDatumIndices['armplane4']], cells=pickedCells)

#Partitioning payload holder flanges
####RMW UPDATED 12-2-2021
p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=BASE_LENGTH+ARM_LENGTH+SPOON_FLANGE_THICKNESS)
armDatumIndices['yz_arm_flange_thickness'] = p.datums.keys()[-1]
p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=BASE_LENGTH+ARM_LENGTH+SPOON_LENGTH-SPOON_FLANGE_THICKNESS)
armDatumIndices['vertical_yz_edgeofspoon'] = p.datums.keys()[-1]
p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=TIP_HEIGHT/2)
armDatumIndices['horizontal_long_spoon_xz_cut'] = p.datums.keys()[-1]
c = p.cells
pickedCells = c.findAt(((BASE_LENGTH+ARM_LENGTH+SPOON_LENGTH/2, 0.0, BASE_WIDTH/2), ))
p.PartitionCellByDatumPlane(datumPlane=d[armDatumIndices['horizontal_long_spoon_xz_cut']], cells=pickedCells)

pickedCells = c.findAt(((BASE_LENGTH+ARM_LENGTH+SPOON_FLANGE_THICKNESS/2, TIP_HEIGHT/2+SPOON_FLANGE_THICKNESS/2, BASE_WIDTH/2), ))
p.PartitionCellByDatumPlane(datumPlane=d[armDatumIndices['yz_arm_flange_thickness']], cells=pickedCells)

pickedCells = c.findAt(((BASE_LENGTH+ARM_LENGTH+SPOON_LENGTH-SPOON_FLANGE_THICKNESS/2, TIP_HEIGHT/2+SPOON_FLANGE_THICKNESS/2, BASE_WIDTH/2), ))
p.PartitionCellByDatumPlane(datumPlane=d[armDatumIndices['vertical_yz_edgeofspoon']], cells=pickedCells)

p = mdb.models['Model-1'].parts['ARM']
p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=0.0)
p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=0.0)
p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=BASE_HEIGHT/2)

#Partitioning spoon pull point node

p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=BASE_LENGTH+ARM_LENGTH+SPOON_LENGTH/2)
armDatumIndices['XZ_spoon_midplane'] = p.datums.keys()[-1]
p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=BASE_WIDTH/2)
armDatumIndices['XY_midplane'] = p.datums.keys()[-1]

pickedCells = c.findAt(((BASE_LENGTH+ARM_LENGTH+SPOON_LENGTH-SPOON_FLANGE_THICKNESS/2, 0.0, BASE_WIDTH/2), ))
p.PartitionCellByDatumPlane(datumPlane=d[armDatumIndices['XY_midplane']], cells=pickedCells)

pickedCells = c.findAt(((BASE_LENGTH+ARM_LENGTH+SPOON_LENGTH-SPOON_FLANGE_THICKNESS/2, 0.0, BASE_WIDTH/4), ),((BASE_LENGTH+ARM_LENGTH+SPOON_LENGTH-SPOON_FLANGE_THICKNESS/2, 0.0, 3*BASE_WIDTH/4), ))
p.PartitionCellByDatumPlane(datumPlane=d[armDatumIndices['XZ_spoon_midplane']], cells=pickedCells)

#Creating Pull Point Set
v = p.vertices
verts = v.findAt(((BASE_LENGTH+ARM_LENGTH+SPOON_LENGTH/2, -TIP_HEIGHT/2, BASE_WIDTH/2), ))
p.Set(vertices=verts, name='SPOON_PULL_POINT')
###RMW CHANGES END

####AXLE#####
#Creating Datum Planes for the AXLE (Connector_beam)
p = mdb.models['Model-1'].parts['Connector_beam']
p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=WALL_THICKNESS/2)
p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=AXLE_LENGTH-WALL_THICKNESS/2)
# p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=0.0)
# p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=AXLE_LENGTH) ##Extraneous planes??
e = p.edges
p = mdb.models['Model-1'].parts['Connector_beam']
p.DatumPointByCoordinate(coords=(0.0, AXLE_DIAMETER/2, AXLE_LENGTH))
p = mdb.models['Model-1'].parts['Connector_beam']
p.DatumPointByCoordinate(coords=(0.0, AXLE_DIAMETER/2, 0.0))
p = mdb.models['Model-1'].parts['Connector_beam']
p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=AXLE_LENGTH/2-BASE_WIDTH/2)
p = mdb.models['Model-1'].parts['Connector_beam']
v = p.vertices
p.DatumPlaneByThreePoints(point1=v.findAt(coordinates=(-AXLE_DIAMETER/4 - OFFSET, AXLE_DIAMETER/4 + OFFSET, 
    AXLE_LENGTH/2 + BASE_WIDTH/2)), point2=v.findAt(coordinates=(-AXLE_DIAMETER/4 - OFFSET, AXLE_DIAMETER/4 + OFFSET, 
    AXLE_LENGTH/2 - BASE_WIDTH/2)), 
    point3=v.findAt(coordinates=(AXLE_DIAMETER/4 - OFFSET, 3*AXLE_DIAMETER/4 + OFFSET, 
    AXLE_LENGTH/2 - BASE_WIDTH/2)))
v2 = p.vertices
p.DatumPlaneByThreePoints(point1=v2.findAt(coordinates=(-AXLE_DIAMETER/2, AXLE_DIAMETER/2, AXLE_LENGTH)), 
    point2=v2.findAt(coordinates=(0.0, 0.0, AXLE_LENGTH)), point3=v2.findAt(
    coordinates=(0.0, 0.0, 0.0)))
v1 = p.vertices
p = mdb.models['Model-1'].parts['Connector_beam']
v = p.vertices
p.DatumPlaneByThreePoints(point1=v.findAt(coordinates=(-AXLE_DIAMETER/2, AXLE_DIAMETER/2, AXLE_LENGTH)), 
    point2=v.findAt(coordinates=(0.0, AXLE_DIAMETER, AXLE_LENGTH)), point3=v.findAt(coordinates=(
    0.0, AXLE_DIAMETER, 0.0)))
p = mdb.models['Model-1'].parts['Connector_beam']
v = p.vertices
p.DatumPlaneByThreePoints(point1=v.findAt(coordinates=(-AXLE_DIAMETER/4 + OFFSET, AXLE_DIAMETER/4 - OFFSET, 
    AXLE_LENGTH/2 + BASE_WIDTH/2)), point2=v.findAt(coordinates=(-AXLE_DIAMETER/4 + OFFSET, AXLE_DIAMETER/4 - OFFSET, 
    AXLE_LENGTH/2 - BASE_WIDTH/2)), 
    point3=v.findAt(coordinates=(AXLE_DIAMETER/4 + OFFSET, 3*AXLE_DIAMETER/4 - OFFSET, 
    AXLE_LENGTH/2 - BASE_WIDTH/2)))
p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=AXLE_LENGTH/2 + BASE_WIDTH/2)

#Partitioning the AXLE
c = p.cells
pickedCells = c.findAt(((0, 0, 0), ))
d = p.datums
p.PartitionCellByDatumPlane(datumPlane=d[6], cells=pickedCells)

c = p.cells
pickedCells = c.findAt(((0, 0, 0), ))
d1 = p.datums
p.PartitionCellByDatumPlane(datumPlane=d1[5], cells=pickedCells)

p = mdb.models['Model-1'].parts['Connector_beam']
c = p.cells
pickedCells = c.findAt(((0, 0, AXLE_LENGTH/2), ))
d = p.datums
p.PartitionCellByDatumPlane(datumPlane=d[9], cells=pickedCells)

c = p.cells
pickedCells = c.findAt(((0, 0, AXLE_LENGTH/2), ))
d1 = p.datums
p.PartitionCellByDatumPlane(datumPlane=d1[14], cells=pickedCells)

c = p.cells
pickedCells = c.findAt(((0, 0, AXLE_LENGTH-0.0001), ), ((0, 
    0, AXLE_LENGTH-WALL_THICKNESS/2-0.0001), ), ((0, 0, AXLE_LENGTH/2-BASE_WIDTH/2-0.0001), ), ((0, 
    0, 0), ))
d = p.datums
p.PartitionCellByDatumPlane(datumPlane=d[10], cells=pickedCells)

c = p.cells
pickedCells = c.findAt(((0, 0, AXLE_LENGTH-0.0001), ), ((0, 
    0, AXLE_LENGTH-WALL_THICKNESS/2-0.0001), ), ((0, 0, AXLE_LENGTH/2-BASE_WIDTH/2-0.0001), ), ((0, 
    0, 0), ))
d1 = p.datums
p.PartitionCellByDatumPlane(datumPlane=d1[13], cells=pickedCells)

#### SIDE WALL 1 ####
#Creating the Datum Points for Side Wall 1
p = mdb.models['Model-1'].parts['Side_wall_1']
p.DatumPointByCoordinate(coords=(WALL_SIDE_X, 0.0, 0.0))
p.DatumPointByCoordinate(coords=(CLEVIS_DIST, 0.0, 0.0))
#Creating the Datum Planes for Side Wall 1
p = mdb.models['Model-1'].parts['Side_wall_1']
p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=0.0)
p = mdb.models['Model-1'].parts['Side_wall_1']
p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=CLEVIS_DIST)
p = mdb.models['Model-1'].parts['Side_wall_1']
p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=WALL_SIDE_X)
p = mdb.models['Model-1'].parts['Side_wall_1']
p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=X1)
p = mdb.models['Model-1'].parts['Side_wall_1']
p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=X2)
p = mdb.models['Model-1'].parts['Side_wall_1']
p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=X4)
p = mdb.models['Model-1'].parts['Side_wall_1']
p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=0.0)
p = mdb.models['Model-1'].parts['Side_wall_1']
p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=P1_Y)
p = mdb.models['Model-1'].parts['Side_wall_1']
p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=P3_Y)
p = mdb.models['Model-1'].parts['Side_wall_1']
p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=-WALL_LENGTH/5)

#Partitioning SIDE WALL 1
c = p.cells
pickedCells = c.findAt(((-WALL_LENGTH/2 + CLEVIS_EDGE_THICK/2, 0.0, WALL_THICKNESS/2), ))
d1 = p.datums
p.PartitionCellByDatumPlane(datumPlane=d1[6], cells=pickedCells)
p = mdb.models['Model-1'].parts['Side_wall_1']
c = p.cells
pickedCells = c.findAt(((WALL_LENGTH/2 - XBEAM_EDGE_THICK/2, 0.0, WALL_THICKNESS/2), ))
d2 = p.datums
p.PartitionCellByDatumPlane(datumPlane=d2[10], cells=pickedCells)
p = mdb.models['Model-1'].parts['Side_wall_1']
c = p.cells
pickedCells = c.findAt((((X4+X3)/2, 0.0, WALL_THICKNESS/2), ))
d1 = p.datums
p.PartitionCellByDatumPlane(datumPlane=d1[9], cells=pickedCells)
p = mdb.models['Model-1'].parts['Side_wall_1']
c = p.cells
pickedCells = c.findAt(((WALL_LENGTH/2 - XBEAM_EDGE_THICK/2, 0.0, WALL_THICKNESS/2), ))
d2 = p.datums
p.PartitionCellByDatumPlane(datumPlane=d2[8], cells=pickedCells)
p = mdb.models['Model-1'].parts['Side_wall_1']
c = p.cells
pickedCells = c.findAt(((WALL_LENGTH/2 - XBEAM_EDGE_THICK/2, 0.0, WALL_THICKNESS/2), ))
d1 = p.datums
p.PartitionCellByDatumPlane(datumPlane=d1[7], cells=pickedCells)
p = mdb.models['Model-1'].parts['Side_wall_1']
c = p.cells
pickedCells = c.findAt(((-WALL_LENGTH/2 + CLEVIS_EDGE_THICK/2, 0.0, WALL_THICKNESS/2), ), ((X4-0.0001, 0.0, WALL_THICKNESS/2), ), (((X4+X3)/2, 0.0, WALL_THICKNESS/2), ), ((WALL_LENGTH/2 - XBEAM_EDGE_THICK/2, 0.0, WALL_THICKNESS/2), ), ((WALL_LENGTH/2 -AXLE_DIAMETER - 3*XBEAM_EDGE_THICK/2, 0.0, WALL_THICKNESS/2), ))
d2 = p.datums
p.PartitionCellByDatumPlane(datumPlane=d2[11], cells=pickedCells)
p = mdb.models['Model-1'].parts['Side_wall_1']
c = p.cells
pickedCells = c.findAt(((-WALL_LENGTH/5, -0.001, WALL_THICKNESS/2), ), ((-WALL_LENGTH/5, 0.001, 
    WALL_THICKNESS/2), ))
d2 = p.datums
p.PartitionCellByDatumPlane(datumPlane=d2[14], cells=pickedCells)
p = mdb.models['Model-1'].parts['Side_wall_1']
c = p.cells
pickedCells = c.findAt(((WALL_LENGTH/2 - XBEAM_EDGE_THICK - AXLE_DIAMETER, -0.001, WALL_THICKNESS/2), ))
d2 = p.datums
p.PartitionCellByDatumPlane(datumPlane=d2[12], cells=pickedCells)
p = mdb.models['Model-1'].parts['Side_wall_1']
c = p.cells
pickedCells = c.findAt(((WALL_LENGTH/2 - XBEAM_EDGE_THICK - AXLE_DIAMETER, 0.001, WALL_THICKNESS/2), ))
d1 = p.datums
p.PartitionCellByDatumPlane(datumPlane=d1[13], cells=pickedCells)
p = mdb.models['Model-1'].parts['Side_wall_1']
v = p.vertices
p.DatumPlaneByThreePoints(point1=v.findAt(coordinates=(P2_X, 0.0, WALL_THICKNESS)), 
    point2=v.findAt(coordinates=(P2_X, 0.0, 0.0)), point3=v.findAt(
    coordinates=(P3_X, P3_Y, 0.0)))
p = mdb.models['Model-1'].parts['Side_wall_1']
v1 = p.vertices
p.DatumPlaneByThreePoints(point1=v1.findAt(coordinates=(P2_X, 0.0, WALL_THICKNESS)), 
    point2=v1.findAt(coordinates=(P2_X, 0.0, 0.0)), point3=v1.findAt(
    coordinates=(P1_X, P1_Y, 0.0)))

# #Creating Pull Point Set
# v = p.vertices
# p.ReferencePoint(point=d1[4])
# mdb.models['Model-1'].parts['Side_wall_1'].features.changeKey(fromName='RP', 
    # toName='AXLE_CENTER_POINT')

#Creating Surfaces on SidelWall 1
#squareNotch
p = mdb.models['Model-1'].parts['Side_wall_1']
s = p.faces
side1Faces = s.findAt(((P1_X+0.01, P1_Y+0.01, WALL_THICKNESS/2), ), ((P1_X-0.01, P1_Y+0.01, WALL_THICKNESS/2), ), ((P3_X-0.01, P3_Y-0.01, WALL_THICKNESS/2), ), ((P3_X+0.01, P3_Y-0.01, WALL_THICKNESS/2), ))
p.Surface(side1Faces=side1Faces, name='SquareNotch')
#circularNotch
s = p.faces
side1Faces = s.findAt(((CLEVIS_DIST + Clev_Opp, Clev_Opp, WALL_THICKNESS/2), ), ((CLEVIS_DIST - Clev_Opp, Clev_Opp, WALL_THICKNESS/2), ), ((CLEVIS_DIST + Clev_Opp, -Clev_Opp, WALL_THICKNESS/2), ), ((CLEVIS_DIST - Clev_Opp, -Clev_Opp, WALL_THICKNESS/2), ))
p.Surface(side1Faces=side1Faces, name='CircularNotch')
#set of the bottom for Boundry Condition
f = p.faces
faces = f.findAt(((X4-0.01, -WALL_HEIGHT/2, WALL_THICKNESS/2), ))
p.Set(faces=faces, name='Bottom')

#### SIDE WALL 2 will now be created ####
p1 = mdb.models['Model-1'].parts['Side_wall_1']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
p = mdb.models['Model-1'].Part(name='Side_wall_2', 
    objectToCopy=mdb.models['Model-1'].parts['Side_wall_1'])
p = mdb.models['Model-1'].parts['Side_wall_2']
p.features['Datum pt-1'].setValues(zValue=WALL_THICKNESS)
p = mdb.models['Model-1'].parts['Side_wall_2']
p.regenerate()
p = mdb.models['Model-1'].parts['Side_wall_2']
p.features['Datum pt-2'].setValues(zValue=WALL_THICKNESS)
p = mdb.models['Model-1'].parts['Side_wall_2']
p.regenerate()



#### CLEVIS ####
#Creating Datum Point For the CLEVIS
p = mdb.models['Model-1'].parts[XBEAM_NAME]
p.DatumPointByCoordinate(coords=((AXLE_LENGTH-WALL_SEP_LENGTH)/2+WALL_SEP_LENGTH/2, XBEAM_HEIGHT, XBEAM_THICKNESS/2))

#Partitioning the CLEVIS
c = p.cells
pickedCells = c.findAt((((AXLE_LENGTH-WALL_SEP_LENGTH)/2+WALL_SEP_LENGTH/2, XBEAM_HEIGHT, XBEAM_THICKNESS/2), ))
e, v1, d1 = p.edges, p.vertices, p.datums
p.PartitionCellByPlanePointNormal(point=d1[6], normal=e.findAt(coordinates=((AXLE_LENGTH-WALL_SEP_LENGTH)/2+WALL_SEP_LENGTH/2-XBEAM_FLAT_LENGTH/2, XBEAM_HEIGHT, XBEAM_THICKNESS/4)), cells=pickedCells)

c = p.cells
pickedCells = c.findAt((((AXLE_LENGTH-WALL_SEP_LENGTH)/2+WALL_SEP_LENGTH/2, XBEAM_HEIGHT, XBEAM_THICKNESS/4), ), (((AXLE_LENGTH-WALL_SEP_LENGTH)/2+WALL_SEP_LENGTH/2, XBEAM_HEIGHT, 3*XBEAM_THICKNESS/4), ))
e1, v, d = p.edges, p.vertices, p.datums
p.PartitionCellByPlanePointNormal(point=d[6], normal=e1.findAt(coordinates=((AXLE_LENGTH-WALL_SEP_LENGTH)/2+WALL_SEP_LENGTH/2+XBEAM_FLAT_LENGTH/4, XBEAM_HEIGHT, 0.0)), cells=pickedCells)
p = mdb.models['Model-1'].parts['CROSSMEMBER']
p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=0.0)
p.DatumPointByCoordinate(coords=(AXLE_LENGTH, XBEAM_AXLE_DIA/2-XBEAM_AXLE_DIA/8, XBEAM_THICKNESS/2))

#Creating Pull Set
v = p.vertices
verts = v.findAt((((AXLE_LENGTH-WALL_SEP_LENGTH)/2+WALL_SEP_LENGTH/2, XBEAM_HEIGHT, XBEAM_THICKNESS/2), ))
p.Set(vertices=verts, name='XBEAM_PULL_POINT')

### New Clevis Partioton for surface grabbing
p = mdb.models['Model-1'].parts['CROSSMEMBER']
p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=WALL_THICKNESS)
p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=AXLE_LENGTH-WALL_THICKNESS)

c = p.cells
pickedCells = c.findAt(((WALL_THICKNESS/2, vertOffset/2, CLEVIS_RAD+XBEAM_THICKNESS/2), ), ((WALL_THICKNESS/2, vertOffset/2, 
    -CLEVIS_RAD+XBEAM_THICKNESS/2), ))
d1 = p.datums
p.PartitionCellByDatumPlane(datumPlane=d1[12], cells=pickedCells)
p = mdb.models['Model-1'].parts['CROSSMEMBER']
c = p.cells
pickedCells = c.findAt(((AXLE_LENGTH-WALL_THICKNESS/2, vertOffset/2, CLEVIS_RAD +XBEAM_THICKNESS/2), ), ((AXLE_LENGTH-WALL_THICKNESS/2, vertOffset/2, 
    -CLEVIS_RAD+XBEAM_THICKNESS/2), ))
d2 = p.datums
p.PartitionCellByDatumPlane(datumPlane=d2[13], cells=pickedCells)

##################################  
#Defining Surfaces for tie constraints (added 11/28/21 by Kirby)
################################## 

#creating surfaces on ARM
a = mdb.models['Model-1'].rootAssembly
p = mdb.models['Model-1'].parts['ARM']
s = p.faces
side1Faces = s.findAt(((0.01, 0, 0), ), ((0.01, BASE_HEIGHT/2, THICKNESS/2), ), ((0.01, BASE_HEIGHT/2, BASE_WIDTH/2), ), ((0.01, BASE_HEIGHT/2, BASE_WIDTH-THICKNESS/2), 
    ), ((0.01, 0, BASE_WIDTH), ), ((0.01, -BASE_HEIGHT/2, THICKNESS/2), ), ((0.01, -BASE_HEIGHT/2, BASE_WIDTH-THICKNESS/2), ), ((0.01, -BASE_HEIGHT/2, BASE_WIDTH/2), ))
p.Surface(side1Faces=side1Faces, name='ArmBase')


#creating surfaces on CROSSMEMBER 
#Left side surface
p = mdb.models['Model-1'].parts['CROSSMEMBER']
s = p.faces
side1Faces = s.findAt(((AXLE_LENGTH-WALL_THICKNESS/2, vertOffset/2, CLEVIS_RAD +XBEAM_THICKNESS/2), ), ((AXLE_LENGTH-WALL_THICKNESS/2, vertOffset/2, 
    -CLEVIS_RAD+XBEAM_THICKNESS/2), ))
p.Surface(side1Faces=side1Faces, name='LeftSide')
#RightSide Surface
s = p.faces
side1Faces = s.findAt(((WALL_THICKNESS/2, vertOffset/2, CLEVIS_RAD+XBEAM_THICKNESS/2), ), ((WALL_THICKNESS/2, vertOffset/2, 
    -CLEVIS_RAD+XBEAM_THICKNESS/2), ))
p.Surface(side1Faces=side1Faces, name='RightSide')
#top of clevis where rope attaches surface
p = mdb.models['Model-1'].parts['CROSSMEMBER']
s = p.faces
zOFF = XBEAM_THICKNESS/4
xOFF = XBEAM_FLAT_LENGTH/4
side1Faces = s.findAt((((AXLE_LENGTH-WALL_SEP_LENGTH)/2+WALL_SEP_LENGTH/2+xOFF, XBEAM_HEIGHT, XBEAM_THICKNESS/2+zOFF), ), (((AXLE_LENGTH-WALL_SEP_LENGTH)/2+WALL_SEP_LENGTH/2-xOFF, XBEAM_HEIGHT, XBEAM_THICKNESS/2+zOFF), ), (((AXLE_LENGTH-WALL_SEP_LENGTH)/2+WALL_SEP_LENGTH/2+xOFF, XBEAM_HEIGHT, XBEAM_THICKNESS/2-zOFF), ), (((AXLE_LENGTH-WALL_SEP_LENGTH)/2+WALL_SEP_LENGTH/2-xOFF, XBEAM_HEIGHT, XBEAM_THICKNESS/2-zOFF), ))
p.Surface(side1Faces=side1Faces, name='CLEVIS_TOP')

#Creating Surfaces on ConnectorBeam
#leftside
p = mdb.models['Model-1'].parts['Connector_beam']
s = p.faces
tempZ = AXLE_LENGTH - 0.0001
AX4 =AXLE_DIAMETER/4
LILOFF = AXLE_DIAMETER/20
AX2 = AXLE_DIAMETER/2
side1Faces = s.findAt(((-AX4, AX2 + AX4, tempZ), ), ((AX4, 
    AX4+AX2, tempZ), ), ((AX4, AX4, tempZ), ), ((-AX4, 
    AX4, tempZ), ), ((-LILOFF, LILOFF, tempZ), ), ((-AX2+LILOFF, 
    AX2 - LILOFF, tempZ), ), ((LILOFF, AXLE_DIAMETER-LILOFF, tempZ), ), ((AX2-LILOFF, 
    AX2+LILOFF, tempZ), ))
p.Surface(side1Faces=side1Faces, name='LeftSide')
#rightSide
s = p.faces
tempZ = 0.0001
side1Faces = s.findAt(((-AX4, AX2 + AX4, tempZ), ), ((AX4, 
    AX4+AX2, tempZ), ), ((AX4, AX4, tempZ), ), ((-AX4, 
    AX4, tempZ), ), ((-LILOFF, LILOFF, tempZ), ), ((-AX2+LILOFF, 
    AX2 - LILOFF, tempZ), ), ((LILOFF, AXLE_DIAMETER-LILOFF, tempZ), ), ((AX2-LILOFF, 
    AX2+LILOFF, tempZ), ))
p.Surface(side1Faces=side1Faces, name='RightSide')
#middleHole
s = p.faces
dumdum = BASE_HEIGHT/2 * cos(pi/4)
side1Faces = s.findAt(((0, AX2, AXLE_LENGTH/2-BASE_WIDTH/2), ), ((0, AX2, AXLE_LENGTH/2+BASE_WIDTH/2), ), ((dumdum, AX2-dumdum, AXLE_LENGTH/2), ), ((-dumdum, AX2+dumdum, AXLE_LENGTH/2), ))
p.Surface(side1Faces=side1Faces, name='MiddleHole')



#Creating Surfaces on SidelWall 1
#squareNotch
p = mdb.models['Model-1'].parts['Side_wall_1']
s = p.faces
side1Faces = s.findAt(((P1_X+0.01, P1_Y+0.01, WALL_THICKNESS/2), ), ((P1_X-0.01, P1_Y+0.01, WALL_THICKNESS/2), ), ((P3_X-0.01, P3_Y-0.01, WALL_THICKNESS/2), ), ((P3_X+0.01, P3_Y-0.01, WALL_THICKNESS/2), ))
p.Surface(side1Faces=side1Faces, name='SquareNotch')
#circularNotch
s = p.faces
side1Faces = s.findAt(((CLEVIS_DIST + Clev_Opp, Clev_Opp, WALL_THICKNESS/2), ), ((CLEVIS_DIST - Clev_Opp, Clev_Opp, WALL_THICKNESS/2), ), ((CLEVIS_DIST + Clev_Opp, -Clev_Opp, WALL_THICKNESS/2), ), ((CLEVIS_DIST - Clev_Opp, -Clev_Opp, WALL_THICKNESS/2), ))
p.Surface(side1Faces=side1Faces, name='CircularNotch')
#set of the bottom for Boundry Condition
f = p.faces
faces = f.findAt(((X4-0.01, -WALL_HEIGHT/2, WALL_THICKNESS/2), ))
p.Set(faces=faces, name='Bottom')

#### SIDE WALL 2 will now be created ####
p1 = mdb.models['Model-1'].parts['Side_wall_1']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
p = mdb.models['Model-1'].Part(name='Side_wall_2', 
    objectToCopy=mdb.models['Model-1'].parts['Side_wall_1'])
p = mdb.models['Model-1'].parts['Side_wall_2']
p.features['Datum pt-1'].setValues(zValue=WALL_THICKNESS)
p = mdb.models['Model-1'].parts['Side_wall_2']
p.regenerate()
p = mdb.models['Model-1'].parts['Side_wall_2']
p.features['Datum pt-2'].setValues(zValue=WALL_THICKNESS)
p = mdb.models['Model-1'].parts['Side_wall_2']
p.regenerate()


##################################  
#Assemble Parts
##################################  
a.regenerate()
a = mdb.models['Model-1'].rootAssembly
del a.features['ARM-1']
a1 = mdb.models['Model-1'].rootAssembly

#Instance bt SIDE WALL 1 and AXLE
a = mdb.models['Model-1'].rootAssembly
p = mdb.models['Model-1'].parts['Side_wall_1']
a.Instance(name='Side_wall_1-1', part=p, dependent=ON)
a = mdb.models['Model-1'].rootAssembly
p = mdb.models['Model-1'].parts['Connector_beam']
a.Instance(name='Connector_beam-1', part=p, dependent=ON)
a1 = mdb.models['Model-1'].rootAssembly
d1 = a1.instances['Side_wall_1-1'].datums
d2 = a1.instances['Connector_beam-1'].datums
a1.FaceToFace(movablePlane=d1[24], fixedPlane=d2[12], flip=OFF, clearance=0.0)
a1 = mdb.models['Model-1'].rootAssembly
d1 = a1.instances['Side_wall_1-1'].datums
d2 = a1.instances['Connector_beam-1'].datums
a1.FaceToFace(movablePlane=d1[25], fixedPlane=d2[11], flip=OFF, clearance=0.0)
a1 = mdb.models['Model-1'].rootAssembly
d1 = a1.instances['Side_wall_1-1'].datums
d2 = a1.instances['Connector_beam-1'].datums
a1.CoincidentPoint(movablePoint=d1[3], fixedPoint=d2[7])

#Instance bt AXLE and SIDE WALL 2
p = mdb.models['Model-1'].parts['Side_wall_2']
a1.Instance(name='Side_wall_2-1', part=p, dependent=ON)
d1 = a1.instances['Side_wall_2-1'].datums
d2 = a1.instances['Connector_beam-1'].datums
a1.FaceToFace(movablePlane=d1[25], fixedPlane=d2[11], flip=OFF, clearance=0.0)
a1 = mdb.models['Model-1'].rootAssembly
d1 = a1.instances['Side_wall_2-1'].datums
d2 = a1.instances['Connector_beam-1'].datums
a1.FaceToFace(movablePlane=d1[24], fixedPlane=d2[12], flip=OFF, clearance=0.0)
a1 = mdb.models['Model-1'].rootAssembly
d1 = a1.instances['Side_wall_2-1'].datums
d2 = a1.instances['Connector_beam-1'].datums
a1.CoincidentPoint(movablePoint=d1[3], fixedPoint=d2[8])



#Instance bt CLEVIS and SIDE WALL 1
p = mdb.models['Model-1'].parts['CROSSMEMBER']
a.Instance(name='CROSSMEMBER-1', part=p, dependent=ON)
a = mdb.models['Model-1'].rootAssembly
f1 = a.instances['CROSSMEMBER-1'].faces
f2 = a.instances['Side_wall_1-1'].faces
a.FaceToFace(movablePlane=f1.findAt(coordinates=(AXLE_LENGTH, 0.106702, 0.060301)), 
    fixedPlane=f2.findAt(coordinates=(1.369204, 0.179423, AXLE_LENGTH)), flip=OFF, 
    clearance=0.0)
a = mdb.models['Model-1'].rootAssembly
d1 = a.instances['CROSSMEMBER-1'].datums
d2 = a.instances['Side_wall_1-1'].datums
a.CoincidentPoint(movablePoint=d1[10], fixedPoint=d2[4])

#Instance bt ARM and AXLE
p = mdb.models['Model-1'].parts['ARM']
a.Instance(name='ARM-1', part=p, dependent=ON)
p1 = a.instances['ARM-1']
p1.translate(vector=(2.334, 0.0, 0.0))
a = mdb.models['Model-1'].rootAssembly
d1 = a.instances['ARM-1'].datums
d2 = a.instances['Side_wall_2-1'].datums
a.FaceToFace(movablePlane=d1[19], fixedPlane=d2[25], flip=OFF, clearance=0.0)
a1 = mdb.models['Model-1'].rootAssembly
d1 = a1.instances['ARM-1'].datums
d2 = a1.instances['Connector_beam-1'].datums
a1.FaceToFace(movablePlane=d1[20], fixedPlane=d2[9], flip=OFF, clearance=0.0)
a1 = mdb.models['Model-1'].rootAssembly
d1 = a1.instances['ARM-1'].datums
d2 = a1.instances['Connector_beam-1'].datums
a1.FaceToFace(movablePlane=d1[21], fixedPlane=d2[10], flip=ON, clearance=0.0)
#: The instance "ARM-1" is fully constrained

# a = mdb.models['Model-1'].rootAssembly
# p = mdb.models['Model-1'].parts['ARM']
# a.Instance(name='ARM-1', part=p, dependent=ON)
# a = mdb.models['Model-1'].rootAssembly
# d1 = a.instances['ARM-1'].datums
# d2 = a.instances['Side_wall_2-1'].datums
# a.FaceToFace(movablePlane=d1[19], fixedPlane=d2[10], flip=ON, clearance=0.0)
# d1 = a.instances['ARM-1'].datums
# d2 = a.instances['Connector_beam-1'].datums
# a.FaceToFace(movablePlane=d1[21], fixedPlane=d2[10], flip=ON, clearance=0.0)
# d1 = a.instances['ARM-1'].datums
# d2 = a.instances['Connector_beam-1'].datums
# a.FaceToFace(movablePlane=d1[20], fixedPlane=d2[9], flip=OFF, clearance=0.0)
# #: The instance "ARM-1" is fully constrained

#extract locations of pull points
# 'AXLE_CENTER_POINT'
spoonCoords = mdb.models['Model-1'].rootAssembly.instances['ARM-1'].sets['SPOON_PULL_POINT'].vertices[0].pointOn[0]
clevisAxleCoords = mdb.models['Model-1'].rootAssembly.instances['Side_wall_1-1'].datums[4].pointOn
print(spoonCoords)
print(clevisAxleCoords)
# spoonCoords = mdb.models['Model-1'].rootAssembly.instances['ARM-1'].sets['SPOON_PULL_POINT'].nodes[0].coordinates

#Rotating the CLEVIS
CLEVIS_ANGLE = 0.0 if abs(spoonCoords[0]-clevisAxleCoords[0]) < 0.001 else atan2(spoonCoords[0]-clevisAxleCoords[0],spoonCoords[1]-clevisAxleCoords[1])*180/pi
print("Updated clevis angle: "+str(CLEVIS_ANGLE))
a = mdb.models['Model-1'].rootAssembly
a1.rotate(instanceList=('CROSSMEMBER-1', ), axisPoint=(WALL_LENGTH - WALL_SIDE_X - CLEVIS_RAD, WALL_HEIGHT/2, WALL_THICKNESS/2), 
    axisDirection=(0.0, 0.0, -AXLE_LENGTH), angle=CLEVIS_ANGLE)




################### Old location of Composite layup




##################################  
#Define Steps
##################################  
print('Defining the Steps')

# Defining Loading step
mdb.models['Model-1'].ImplicitDynamicsStep(name='Loading', previous='Initial', 
    maxNumInc=30, application=QUASI_STATIC, initialInc=0.1, nohaf=OFF, 
    amplitude=RAMP, alpha=DEFAULT, initialConditions=OFF, nlgeom=ON)
# Defining Launch step
mdb.models['Model-1'].ImplicitDynamicsStep(name='Launch', previous='Loading', 
    timePeriod=2.0, application=TRANSIENT_FIDELITY, initialInc=0.2, 
    minInc=2e-05, nohaf=OFF, initialConditions=ON)



##################################  
#Create Interactions
##################################  
print('Defining Interactions')

#Tieing Arm to the connector beam
region1=a.instances['ARM-1'].surfaces['ArmBase']
region2=a.instances['Connector_beam-1'].surfaces['MiddleHole']
mdb.models['Model-1'].Tie(name='Arm__Beam', master=region1, slave=region2, 
    positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)

#Tieing Connector Beam to SideWalls
region1=a.instances['Connector_beam-1'].surfaces['LeftSide']
region2=a.instances['Side_wall_1-1'].surfaces['SquareNotch']
mdb.models['Model-1'].Tie(name='Beam__Sidewall1', master=region1, slave=region2, 
    positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)
region1=a.instances['Connector_beam-1'].surfaces['RightSide']
a = mdb.models['Model-1'].rootAssembly
region2=a.instances['Side_wall_2-1'].surfaces['SquareNotch']
mdb.models['Model-1'].Tie(name='Beam__Sidewall2', master=region1, 
    slave=region2, positionToleranceMethod=COMPUTED, adjust=ON, 
    tieRotations=ON, thickness=ON)
	
#Tieing CrossBeam to sidewalls
region1=a.instances['CROSSMEMBER-1'].surfaces['LeftSide']
region2=a.instances['Side_wall_1-1'].surfaces['CircularNotch']
mdb.models['Model-1'].Tie(name='Cross__Sidewall1', master=region1, 
    slave=region2, positionToleranceMethod=COMPUTED, adjust=ON, 
    tieRotations=ON, thickness=ON)
region1=a.instances['CROSSMEMBER-1'].surfaces['RightSide']
region2=a.instances['Side_wall_2-1'].surfaces['CircularNotch']
mdb.models['Model-1'].Tie(name='Cross__Sidewall2', master=region1, 
    slave=region2, positionToleranceMethod=COMPUTED, adjust=ON, 
    tieRotations=ON, thickness=ON)
    
#Creating Torsional Spring connector 
mdb.models['Model-1'].ConnectorSection(name='ConnSect-1', 
    rotationalType=ROTATION)
elastic_0 = connectorBehavior.ConnectorElasticity(components=(4, ), table=((
    420.0, ), )) #Change this value for stiffness coefficient 
mdb.models['Model-1'].sections['ConnSect-1'].setValues(behaviorOptions =(
    elastic_0, ) )
mdb.models['Model-1'].sections['ConnSect-1'].behaviorOptions[0].ConnectorOptions(
    )
a = mdb.models['Model-1'].rootAssembly
d1 = a.instances['Side_wall_2-1'].datums
rp1 = a.ReferencePoint(point=d1[7])
rpid1 = a.referencePoints[rp1.id]
oldName = rp1.name
mdb.models['Model-1'].rootAssembly.features.changeKey(fromName=oldName, 
    toName='RP-1')
e1 = a.instances['Side_wall_1-1'].edges
rp2 = a.ReferencePoint(point=a.instances['Side_wall_1-1'].InterestingPoint(
    edge=e1.findAt(coordinates=(0.475, 0.2, 0.05)), rule=MIDDLE))
rpid2 = a.referencePoints[rp2.id]
oldName = rp2.name
mdb.models['Model-1'].rootAssembly.features.changeKey(fromName=oldName, 
    toName='RP-2')
dtm1 = a.DatumCsysByThreePoints(origin=d1[7], coordSysType=CARTESIAN, 
    point1=a.instances['Side_wall_1-1'].InterestingPoint(edge=e1.findAt(
    coordinates=(0.475, 0.2, 0.05)), rule=MIDDLE), isZ=True)
dtmid1 = a.datums[dtm1.id]
a = mdb.models['Model-1'].rootAssembly
wire = a.WirePolyLine(points=((rpid1, rpid2), ), mergeType=IMPRINT, 
    meshable=False)
oldName = wire.name
mdb.models['Model-1'].rootAssembly.features.changeKey(fromName=oldName, 
    toName='Wire-1')
a = mdb.models['Model-1'].rootAssembly
e1 = a.edges
edges1 = e1.findAt(((0.4875, 0.2, -0.25), ))
a.Set(edges=edges1, name='Wire-1-Set-1')
region = mdb.models['Model-1'].rootAssembly.sets['Wire-1-Set-1']
csa = a.SectionAssignment(sectionName='ConnSect-1', region=region)
#: The section "ConnSect-1" has been assigned to 1 wire or attachment line.
a.ConnectorOrientation(region=csa.getSet(), localCsys1=dtmid1)



##################################  
#Create Loads
##################################  
print('Defining Loads')

    
    
    
    
##################################  
#Define BCs
##################################  
print('Defining all BCs')
#Clamping the bottom of LunaCat to the floor
region = a.instances['Side_wall_1-1'].sets['Bottom']
mdb.models['Model-1'].EncastreBC(name='Sidewall1Clamp', 
    createStepName='Initial', region=region, localCsys=None)
region = a.instances['Side_wall_2-1'].sets['Bottom']
mdb.models['Model-1'].EncastreBC(name='Sidewall2Clamp', 
    createStepName='Initial', region=region, localCsys=None)




##################################  
#Define Sets
##################################  
print('Defining Sets')
########### Start of Temp Code #################
# create set of bottom of arm
p = mdb.models['Model-1'].parts['ARM']
f = p.faces
faces = f.findAt(((ARM_LENGTH+WALL_SIDE_X, 0.0, 0.0), ))
p.Set(faces=faces, name='TipOfArm')

# displace bottom of arm for loading and release for launch
region = a.instances['ARM-1'].sets['TipOfArm']
mdb.models['Model-1'].DisplacementBC(name='Winch', createStepName='Loading', 
    region=region, u1=0.1, u2=-0.5, u3=0.0, ur1=0.0, ur2=0.0, ur3=UNSET, 
    amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
    localCsys=None)
mdb.models['Model-1'].boundaryConditions['Winch'].deactivate('Launch')



##################################      
#Mesh Parts
##################################  
print('Meshing the Part')
###Temporary meshing ARM
p = mdb.models['Model-1'].parts['ARM']
c = p.cells
pickedCells = c.findAt(((2.395333, 0.012, 0.046667), ))
d = p.datums
p.PartitionCellByDatumPlane(datumPlane=d[15], cells=pickedCells)
p.seedPart(size=0.05, deviationFactor=0.1, minSizeFactor=0.1)
p = mdb.models['Model-1'].parts['ARM']
p.generateMesh()
###Temporary meshing ConnectorBeam
# p = mdb.models['Model-1'].parts['Connector_beam']
# p.seedPart(size=0.044, deviationFactor=0.1, minSizeFactor=0.1)
# p.generateMesh()
###Temporary meshing CROSSBEAM
p = mdb.models['Model-1'].parts['CROSSMEMBER']
c = p.cells
pickedRegions = c.findAt(((2.0, 0.105796, 0.037244), ), ((0.0, 0.105796, 
    0.037244), ), ((0.0, 0.105796, 0.062756), ), ((2.0, 0.105796, 0.062756), ))
p.setMeshControls(regions=pickedRegions, elemShape=TET, technique=FREE)
elemType1 = mesh.ElemType(elemCode=C3D20R)
elemType2 = mesh.ElemType(elemCode=C3D15)
elemType3 = mesh.ElemType(elemCode=C3D10)
p = mdb.models['Model-1'].parts['CROSSMEMBER']
c = p.cells
cells = c.findAt(((2.0, 0.105796, 0.037244), ), ((0.0, 0.105796, 0.037244), ), 
    ((0.0, 0.105796, 0.062756), ), ((2.0, 0.105796, 0.062756), ))
pickedRegions =(cells, )
p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, 
    elemType3))
p = mdb.models['Model-1'].parts['CROSSMEMBER']
p.seedPart(size=0.052, deviationFactor=0.1, minSizeFactor=0.1)
p = mdb.models['Model-1'].parts['CROSSMEMBER']
p.generateMesh()
###Temporary meshing SideWalls
p = mdb.models['Model-1'].parts['Side_wall_1']
c = p.cells
pickedRegions = c.findAt(((0.433333, 0.166667, 0.033333), ), ((0.466667, 
    0.166667, 0.05), ), ((0.466667, 0.266667, 0.033333), ), ((0.666667, 0.0, 
    0.016667), ), ((0.533333, 0.166667, 0.05), ), ((0.666667, 0.4, 0.033333), 
    ), ((0.466667, 0.233333, 0.05), ), ((1.487022, 0.100846, 0.016667), ), ((
    1.487244, 0.135871, 0.05), ), ((1.512978, 0.100846, 0.033333), ), ((
    1.933333, 0.066667, 0.066667), ), ((1.512978, 0.299154, 0.016667), ), ((
    1.333333, 0.333333, 0.1), ), ((0.533333, 0.233333, 0.05), ), ((1.512756, 
    0.264129, 0.05), ))
p.setMeshControls(regions=pickedRegions, elemShape=TET, technique=FREE)
elemType1 = mesh.ElemType(elemCode=C3D20R)
elemType2 = mesh.ElemType(elemCode=C3D15)
elemType3 = mesh.ElemType(elemCode=C3D10)
p = mdb.models['Model-1'].parts['Side_wall_1']
c = p.cells
cells = c.findAt(((0.433333, 0.166667, 0.033333), ), ((0.466667, 0.166667, 
    0.05), ), ((0.466667, 0.266667, 0.033333), ), ((0.666667, 0.0, 0.016667), 
    ), ((0.533333, 0.166667, 0.05), ), ((0.666667, 0.4, 0.033333), ), ((
    0.466667, 0.233333, 0.05), ), ((1.487022, 0.100846, 0.016667), ), ((
    1.487244, 0.135871, 0.05), ), ((1.512978, 0.100846, 0.033333), ), ((
    1.933333, 0.066667, 0.066667), ), ((1.512978, 0.299154, 0.016667), ), ((
    1.333333, 0.333333, 0.1), ), ((0.533333, 0.233333, 0.05), ), ((1.512756, 
    0.264129, 0.05), ))
pickedRegions =(cells, )
p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, 
    elemType3))
p = mdb.models['Model-1'].parts['Side_wall_1']
c = p.cells
pickedRegions = c.findAt(((1.487022, 0.299154, 0.033333), ))
p.setMeshControls(regions=pickedRegions, elemShape=TET, technique=FREE)
elemType1 = mesh.ElemType(elemCode=C3D20R)
elemType2 = mesh.ElemType(elemCode=C3D15)
elemType3 = mesh.ElemType(elemCode=C3D10)
p = mdb.models['Model-1'].parts['Side_wall_1']
c = p.cells
cells = c.findAt(((1.487022, 0.299154, 0.033333), ))
pickedRegions =(cells, )
p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, 
    elemType3))
p = mdb.models['Model-1'].parts['Side_wall_1']
p.seedPart(size=0.05, deviationFactor=0.1, minSizeFactor=0.1)
p = mdb.models['Model-1'].parts['Side_wall_1']
p.generateMesh()
p = mdb.models['Model-1'].parts['Side_wall_2']
c = p.cells
pickedRegions = c.findAt(((0.466667, 0.133333, 0.066667), ), ((0.466667, 
    0.166667, 0.05), ), ((0.433333, 0.233333, 0.066667), ), ((0.666667, 0.0, 
    0.066667), ), ((0.533333, 0.166667, 0.05), ), ((0.533333, 0.266667, 
    0.066667), ), ((0.666667, 0.4, 0.033333), ), ((1.487022, 0.100846, 
    0.066667), ), ((1.487244, 0.135871, 0.05), ), ((1.333333, 0.4, 0.066667), 
    ), ((1.512978, 0.100846, 0.083333), ), ((1.866667, 0.133333, 0.033333), ), 
    ((1.512978, 0.299154, 0.066667), ), ((1.487244, 0.264129, 0.05), ), ((
    0.466667, 0.4, 0.016667), ), ((1.733333, 0.266667, 0.016667), ))
p.setMeshControls(regions=pickedRegions, elemShape=TET, technique=FREE)
elemType1 = mesh.ElemType(elemCode=C3D20R)
elemType2 = mesh.ElemType(elemCode=C3D15)
elemType3 = mesh.ElemType(elemCode=C3D10)
c = p.cells
cells = c.findAt(((0.466667, 0.133333, 0.066667), ), ((0.466667, 0.166667, 
    0.05), ), ((0.433333, 0.233333, 0.066667), ), ((0.666667, 0.0, 0.066667), 
    ), ((0.533333, 0.166667, 0.05), ), ((0.533333, 0.266667, 0.066667), ), ((
    0.666667, 0.4, 0.033333), ), ((1.487022, 0.100846, 0.066667), ), ((
    1.487244, 0.135871, 0.05), ), ((1.333333, 0.4, 0.066667), ), ((1.512978, 
    0.100846, 0.083333), ), ((1.866667, 0.133333, 0.033333), ), ((1.512978, 
    0.299154, 0.066667), ), ((1.487244, 0.264129, 0.05), ), ((0.466667, 0.4, 
    0.016667), ), ((1.733333, 0.266667, 0.016667), ))
pickedRegions =(cells, )
p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, 
    elemType3))
p.seedPart(size=0.05, deviationFactor=0.1, minSizeFactor=0.1)
p.generateMesh()
p = mdb.models['Model-1'].parts['ARM']
del mdb.models['Model-1'].parts['ARM'].sectionAssignments[0]
p = mdb.models['Model-1'].parts['ARM']
region = p.sets['Set-1']
p = mdb.models['Model-1'].parts['ARM']
p.SectionAssignment(region=region, sectionName='Section-1', offset=0.0, 
    offsetType=MIDDLE_SURFACE, offsetField='', 
    thicknessAssignment=FROM_SECTION)
del mdb.models['Model-1'].materials['Graphite Epoxy AS/3501']
del mdb.models['Model-1'].sections['Section-2']

#meshing the axle
p = mdb.models['Model-1'].parts['Connector_beam']
p.seedPart(size=0.014, deviationFactor=0.1, minSizeFactor=0.1)
p = mdb.models['Model-1'].parts['Connector_beam']
p.generateMesh()
STTAAHHHPPPP



#here

    # ####################################
    # ## Creation/Execution of the Job ###
    # ####################################
    # print 'Creating/Running Job'

    # ModelName='Model-1'

    # mdb.Job(name=ModelName, model=ModelName, description='', type=ANALYSIS, 
        # atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
        # memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
        # explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
        # modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
        # scratch='', multiprocessingMode=DEFAULT, numCpus=1, numGPUs=0)

    # job=mdb.jobs[ModelName]

    # delete lock file, which for some reason tends to hang around, if it exists
    # if os.access('%s.lck'%ModelName,os.F_OK):
        # os.remove('%s.lck'%ModelName)
        
    # Run the job, then process the results.        
    # job.submit()
    # job.waitForCompletion()
    # print 'Completed job'

    # tipDisp = getResults(ModelName)
    # DataFile = open('PostData.txt','a')
    # DataFile.write('%10f %10f\n' % (sizes,tipDisp))
    # DataFile.close()	


    # Calculations (if needed)


    # print sizes
    # ##END LOOP (i.e., end indentation)

print('DONE!!')
