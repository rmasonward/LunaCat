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
import Post_P_Script_Velo
import Post_P_Script

def evalModel(BASE_WIDTH, BASE_HEIGHT, THICKNESS, ARM_LENGTH, TAPER_RATIO, WALL_LENGTH, AXLE_LENGTH, MATERIAL_TYPE, ToptThickness, L1_percent, L2_percent, CLEVIS_EDGE_THICK, XBEAM_EDGE_THICK):
    print(BASE_WIDTH, BASE_HEIGHT, THICKNESS, ARM_LENGTH, TAPER_RATIO, WALL_LENGTH, AXLE_LENGTH, MATERIAL_TYPE, ToptThickness, L1_percent, L2_percent, CLEVIS_EDGE_THICK, XBEAM_EDGE_THICK)
    ######################################
    # Variable and Fixed Design Parameters
    ######################################

    #Arm
    #BASE_WIDTH = 0.15 #m ###DESIGN VARIABLE
    #BASE_HEIGHT = 0.07 #m ###DESIGN VARIABLE
    #TAPER_RATIO = 3.5 #DESIGN VARIABLE
    BASE_LENGTH = 0.2 #m
    TIP_HEIGHT = TAPER_RATIO * BASE_HEIGHT #m


    #THICKNESS = 0.004 #m ###DESIGN VARIABLE
    if THICKNESS*2 > min(BASE_HEIGHT,BASE_WIDTH,TIP_HEIGHT):
        raise ValueError("Thickness is too big")

    #ARM_LENGTH = 3.7 #m ###DESIGN VARIABLE
    SPOON_LENGTH = 0.2 #m
    SPOON_FLANGE_THICKNESS = 0.002 #m
    SPOON_FLANGE_HEIGHT = 0.01 #m

    ## PAYLOAD VARIABLES
    PAY_DESNITY = 9000

    WALL_THICKNESS = 0.1#m
    #WALL_LENGTH = 3.0 #m #DESIGN VARIABLE
    WALL_HEIGHT = 0.4 #m

    DRAW_WIRE_RADIUS = 0.01
    WIRE_DRAW_DISTANCE = 0.2

    # Axle Beam Variables
    S1 = 0.1 #m  used to sketch the cross section
    #AXLE_LENGTH = 1.2 #m ###DESIGN VARIABLE
    AXLE_DIAMETER = 0.2 #m  this is defined as the hypotenuse of the square cross section


    # Cross Member Variables
    WALL_SEP_LENGTH = AXLE_LENGTH-2*WALL_THICKNESS #m

    XBEAM_HEIGHT = 0.3 #m
    XBEAM_FLAT_LENGTH = 0.02 #m
    XBEAM_THICKNESS = 0.1 #m
    XBEAM_AXLE_DIA = 0.2 #m
    OFFSET = sqrt((BASE_HEIGHT**2)/8) 
    CLEVIS_RAD = XBEAM_AXLE_DIA/2
    Clev_Opp = CLEVIS_RAD * cos(pi/4)
    
    
    #Material properties
    if MATERIAL_TYPE == 1.0:
        material_name = 'Aluminum-2024'
        density = 2780.0 #kg/m^3
        youngs = 73100000000.0 #Pa
        poisson = 0.33 
    elif MATERIAL_TYPE == 2.0:
        material_name = 'Aluminum-6061'
        density = 2700.0 #kg/m^3
        youngs = 68900000000.0 #Pa
        poisson = 0.33 
    elif MATERIAL_TYPE == 3.0:
        material_name = 'Ti-6Al-4V'
        density = 4429.0 #kg/m^3
        youngs = 104000000000.0 #Pa
        poisson = 0.31


    # New variables from top opt
    #### 5 design variables for top opt
    #ToptThickness = 0.075 #KDV 1
    L1 = L1_percent/100 * WALL_LENGTH
    L2 = L2_percent/100 * WALL_LENGTH
    #CLEVIS_EDGE_THICK = 0.05 #KDV 4
    # XBEAM_EDGE_THICK = 0.05 #KDV 5

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

    # #Make the payload holders here
    # t = p.MakeSketchTransform(sketchPlane=f.findAt(coordinates=(ARM_LENGTH+BASE_LENGTH-SPOON_LENGTH/2, 0.0, 
        # BASE_WIDTH)), sketchUpEdge=e.findAt(coordinates=(ARM_LENGTH+BASE_LENGTH+SPOON_LENGTH, 0.0, BASE_WIDTH)), 
        # sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(ARM_LENGTH+BASE_LENGTH, TIP_HEIGHT/2, 
        # BASE_WIDTH))
    # s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=11.2, 
        # gridSpacing=0.28, transform=t)
    # g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    # s.sketchOptions.setValues(decimalPlaces=3)
    # s.setPrimaryObject(option=SUPERIMPOSE)
    # p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)

    # #proximal flange
    # s.Line(point1=(0.0,0.0), point2=(0.0, SPOON_FLANGE_HEIGHT))
    # s.Line(point1=(0.0, SPOON_FLANGE_HEIGHT), point2=(SPOON_FLANGE_THICKNESS, SPOON_FLANGE_HEIGHT))
    # s.Line(point1=(SPOON_FLANGE_THICKNESS, SPOON_FLANGE_HEIGHT), point2=(SPOON_FLANGE_THICKNESS, SPOON_FLANGE_THICKNESS))
    # s.Line(point1=(SPOON_FLANGE_THICKNESS, SPOON_FLANGE_THICKNESS), point2=(SPOON_FLANGE_HEIGHT, SPOON_FLANGE_THICKNESS))
    # s.Line(point1=(SPOON_FLANGE_HEIGHT, SPOON_FLANGE_THICKNESS), point2=(SPOON_FLANGE_HEIGHT, 0.0))
    # s.Line(point1=(SPOON_FLANGE_HEIGHT, 0.0), point2=(0.0, 0.0))

    # #distal flange
    # s.Line(point1=(SPOON_LENGTH-SPOON_FLANGE_HEIGHT,0.0), point2=(SPOON_LENGTH-SPOON_FLANGE_HEIGHT,SPOON_FLANGE_THICKNESS))
    # s.Line(point1=(SPOON_LENGTH-SPOON_FLANGE_HEIGHT,SPOON_FLANGE_THICKNESS), point2=(SPOON_LENGTH-SPOON_FLANGE_THICKNESS,SPOON_FLANGE_THICKNESS))
    # s.Line(point1=(SPOON_LENGTH-SPOON_FLANGE_THICKNESS,SPOON_FLANGE_THICKNESS), point2=(SPOON_LENGTH-SPOON_FLANGE_THICKNESS,SPOON_FLANGE_HEIGHT))
    # s.Line(point1=(SPOON_LENGTH-SPOON_FLANGE_THICKNESS,SPOON_FLANGE_HEIGHT), point2=(SPOON_LENGTH,SPOON_FLANGE_HEIGHT))
    # s.Line(point1=(SPOON_LENGTH,SPOON_FLANGE_HEIGHT), point2=(SPOON_LENGTH,0.0))
    # s.Line(point1=(SPOON_LENGTH,0.0), point2=(SPOON_LENGTH-SPOON_FLANGE_HEIGHT,0.0))

    # p.SolidExtrude(sketchPlane=f1.findAt(coordinates=(ARM_LENGTH+BASE_LENGTH-SPOON_LENGTH/2, 0.0, BASE_WIDTH)), 
        # sketchUpEdge=e1.findAt(coordinates=(ARM_LENGTH+BASE_LENGTH+SPOON_LENGTH, 0.0, BASE_WIDTH)), 
        # sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, sketch=s, depth=BASE_WIDTH, 
        # flipExtrudeDirection=ON)
    # s.unsetPrimaryObject()
    # del mdb.models['Model-1'].sketches['__profile__']
    
    # #Deleting flanges
    # p = mdb.models['Model-1'].parts['ARM']
    # del p.features['Solid extrude-2']
    ########################################################################

    ####Making Parameterized payload
    PAY_LENGTH = SPOON_LENGTH - 5.0 * SPOON_FLANGE_THICKNESS
    PAY_HEIGHT = PAY_LENGTH/2
    PAY_WIDTH = BASE_WIDTH
    s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=10.0)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=STANDALONE)
    s.rectangle(point1=(0.0, 0.0), point2=(PAY_LENGTH, PAY_HEIGHT))
    p = mdb.models['Model-1'].Part(name='Payload', dimensionality=THREE_D, 
        type=DEFORMABLE_BODY)
    p = mdb.models['Model-1'].parts['Payload']
    p.BaseSolidExtrude(sketch=s, depth=PAY_WIDTH)
    s.unsetPrimaryObject()
    p = mdb.models['Model-1'].parts['Payload']
    del mdb.models['Model-1'].sketches['__profile__']

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

    #Datum plane on top side of connector beam close to clevis
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

    #Creating Material
    print('Creating the Materials')
    mdb.models['Model-1'].Material(name=material_name)
    mdb.models['Model-1'].materials[material_name].Density(table=((density, ), ))
    mdb.models['Model-1'].materials[material_name].Elastic(table=((youngs, 
        poisson), ))

    # #Creating Graphite
    # mdb.models['Model-1'].Material(name='Graphite Epoxy AS/3501')
    # mdb.models['Model-1'].materials['Graphite Epoxy AS/3501'].Elastic(type=ENGINEERING_CONSTANTS, 
        # table=((137894900000.0, 8963168000.0, 
        # 8963168000.0, 2068423000.0, 2068423000.0, 2068423000.0, 6894745000.0, 
        # 6894745000.0, 6894745000.0), ))
    # mdb.models['Model-1'].materials['Graphite Epoxy AS/3501'].Density(table=((
        # 1610.0, ), ))

        
    #Creating Payload material
    mdb.models['Model-1'].Material(name='MoonRock')
    mdb.models['Model-1'].materials['MoonRock'].Density(table=((PAY_DESNITY, ), ))
    mdb.models['Model-1'].materials['MoonRock'].Elastic(table=((10000000000.0, 
        0.3), ))

    ##################################  
    #Create/Assign Section
    ##################################  

    print('Creating the Sections')
    #Creating section 1
    p = mdb.models['Model-1'].parts['Side_wall_1']
    mdb.models['Model-1'].HomogeneousSolidSection(name='Section-1', 
        material=material_name, thickness=None)
        
    #Creating Section 2
    mdb.models['Model-1'].HomogeneousSolidSection(name='Section-2', 
        material='Graphite Epoxy AS/3501', thickness=None)

    #Creating Section 3
    mdb.models['Model-1'].HomogeneousSolidSection(name='Section-3', 
        material='MoonRock', thickness=None)

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

    #Assigning Section 1 to ARM
    p = mdb.models['Model-1'].parts['ARM']
    c = p.cells
    cells = c.findAt(((0.0, 0.0, 0.0), ))
    region = p.Set(cells=cells, name='Set-1')
    p = mdb.models['Model-1'].parts['ARM']
    p.SectionAssignment(region=region, sectionName='Section-1', offset=0.0, 
        offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)

    #Assigning Section 3 to PAYLOAD
    p = mdb.models['Model-1'].parts['Payload']
    c = p.cells
    cells = c.findAt(((PAY_LENGTH/2, PAY_HEIGHT/2, PAY_WIDTH/2), ))
    region = p.Set(cells=cells, name='Set-1')
    p = mdb.models['Model-1'].parts['Payload']
    p.SectionAssignment(region=region, sectionName='Section-3', offset=0.0, 
        offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)
        
        
    ###DATUM ARRAY HEADQUARTERS
    armDatumIndices = {} ### Creates array to reference datum plane
    connectorBeamDatumIndices = {} ##Array for crossmember
    payloadDatumIndices = {} ##Array for payload
    
    
    ##################################  
    # Create Composite Layup
    ##################################  

    # p = mdb.models['Model-1'].parts['ARM']
    # v, e = p.vertices, p.edges
    # d = p.datums
        
    # p.DatumCsysByThreePoints(point2=v.findAt(coordinates=(0.0, BASE_HEIGHT/2, 0.0)), 
        # name='Datum csys-1', coordSysType=CARTESIAN, origin=(0.0, 0.0, 0.0), 
        # point1=p.InterestingPoint(edge=e.findAt(coordinates=(ARM_LENGTH+BASE_LENGTH+SPOON_LENGTH, 0.0, 0.0)), 
        # rule=MIDDLE))
        
    # armDatumIndices['Composite-csys'] = d.keys()[-1]
    # layupOrientation = d[armDatumIndices['Composite-csys']]

    # c = p.cells
    # cells = c.findAt(((0.0, 0.0, 0.0), ), ((BASE_LENGTH+ARM_LENGTH+SPOON_LENGTH, 0.0, 0.0), 
        # ))
    # region1=regionToolset.Region(cells=cells)
    # region2=regionToolset.Region(cells=cells)
    # region3=regionToolset.Region(cells=cells)
    # compositeLayup = mdb.models['Model-1'].parts['ARM'].CompositeLayup(
        # name='CompositeLayup-1', description='', elementType=SOLID, 
        # symmetric=False, thicknessAssignment=FROM_SECTION)
    # compositeLayup.ReferenceOrientation(orientationType=SYSTEM, 
        # localCsys=layupOrientation, fieldName='', 
        # additionalRotationType=ROTATION_NONE, angle=0.0, 
        # additionalRotationField='', axis=AXIS_3, stackDirection=STACK_3)
    # compositeLayup.CompositePly(suppressed=True, plyName='Ply-1', region=region1, 
        # material='Graphite Epoxy AS/3501', thicknessType=SPECIFY_THICKNESS, 
        # thickness=0.1, orientationType=ANGLE_0, 
        # additionalRotationType=ROTATION_NONE, additionalRotationField='', 
        # axis=AXIS_3, angle=0.0, numIntPoints=3)
    # compositeLayup.CompositePly(suppressed=True, plyName='Ply-2', region=region2, 
        # material='Graphite Epoxy AS/3501', thicknessType=SPECIFY_THICKNESS, 
        # thickness=0.1, orientationType=ANGLE_0, 
        # additionalRotationType=ROTATION_NONE, additionalRotationField='', 
        # axis=AXIS_3, angle=45, numIntPoints=3)
    # compositeLayup.CompositePly(suppressed=True, plyName='Ply-3', region=region3, 
        # material='Graphite Epoxy AS/3501', thicknessType=SPECIFY_THICKNESS, 
        # thickness=0.1, orientationType=ANGLE_0, 
        # additionalRotationType=ROTATION_NONE, additionalRotationField='', 
        # axis=AXIS_3, angle=-45, numIntPoints=3)



    ##################################  
    #Defining the face partitions
    ##################################  
    print('Partitioning part')
    print('Partitioning arm')
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
    armDatumIndices['spoon_terminus_yz'] = p.datums.keys()[-1]
    p.PartitionCellByDatumPlane(datumPlane=d[armDatumIndices['spoon_terminus_yz']], cells=pickedCells)

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
    # pickedCells = c.findAt(((BASE_LENGTH+ARM_LENGTH+SPOON_LENGTH/2, 0.0, BASE_WIDTH/2), ))
    #p.PartitionCellByDatumPlane(datumPlane=d[armDatumIndices['horizontal_long_spoon_xz_cut']], cells=pickedCells)

    # pickedCells = c.findAt(((BASE_LENGTH+ARM_LENGTH+SPOON_FLANGE_THICKNESS/2, TIP_HEIGHT/2+SPOON_FLANGE_THICKNESS/2, BASE_WIDTH/2), ))
    #p.PartitionCellByDatumPlane(datumPlane=d[armDatumIndices['yz_arm_flange_thickness']], cells=pickedCells)

    # pickedCells = c.findAt(((BASE_LENGTH+ARM_LENGTH+SPOON_LENGTH-SPOON_FLANGE_THICKNESS/2, TIP_HEIGHT/2+SPOON_FLANGE_THICKNESS/2, BASE_WIDTH/2), ))
    #p.PartitionCellByDatumPlane(datumPlane=d[armDatumIndices['vertical_yz_edgeofspoon']], cells=pickedCells)
    
    
    #Datum planes for Assembly
    p = mdb.models['Model-1'].parts['ARM']
    p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=0.0)
    armDatumIndices['Arm_base_assembly_plane_yz_a1'] = p.datums.keys()[-1]
    p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=0.0)
    armDatumIndices['Arm_side_principal_plane_xy'] = p.datums.keys()[-1]
    p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=BASE_WIDTH)
    armDatumIndices['Arm_side_assembly_plane_xy_a2'] = p.datums.keys()[-1]
    p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=BASE_HEIGHT/2)
    armDatumIndices['Arm_top_assembly_plane_xz_a3'] = p.datums.keys()[-1]

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
    print('Partitioning axle')
    ####AXLE#####
    #Creating Datum Planes for the AXLE (Connector_beam)
    p = mdb.models['Model-1'].parts['Connector_beam']
    d = p.datums
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
    #Datum plane on top inside side of connector beam slot parallel to arm
    p.DatumPlaneByThreePoints(point1=v.findAt(coordinates=(-AXLE_DIAMETER/4 - OFFSET, AXLE_DIAMETER/4 + OFFSET, 
        AXLE_LENGTH/2 + BASE_WIDTH/2)), point2=v.findAt(coordinates=(-AXLE_DIAMETER/4 - OFFSET, AXLE_DIAMETER/4 + OFFSET, 
        AXLE_LENGTH/2 - BASE_WIDTH/2)), 
        point3=v.findAt(coordinates=(AXLE_DIAMETER/4 - OFFSET, 3*AXLE_DIAMETER/4 + OFFSET, 
        AXLE_LENGTH/2 - BASE_WIDTH/2)))
    connectorBeamDatumIndices['Crossbeam_top_inside_surface_assembly_a3'] = d.keys()[-1]
    v2 = p.vertices
    #Datum plane on bottom surface of connector beam facing the ground orthogonal to arm
    p.DatumPlaneByThreePoints(point1=v2.findAt(coordinates=(-AXLE_DIAMETER/2, AXLE_DIAMETER/2, AXLE_LENGTH)), 
        point2=v2.findAt(coordinates=(0.0, 0.0, AXLE_LENGTH)), point3=v2.findAt(
        coordinates=(0.0, 0.0, 0.0)))
    connectorBeamDatumIndices['Crossbeam_base_surface_assembly_a1'] = d.keys()[-1]
    v1 = p.vertices
    p = mdb.models['Model-1'].parts['Connector_beam']
    v = p.vertices
    #Datum plane on side face of connector beam facing away from catapult parallel to arm
    p.DatumPlaneByThreePoints(point1=v.findAt(coordinates=(-AXLE_DIAMETER/2, AXLE_DIAMETER/2, AXLE_LENGTH)), 
        point2=v.findAt(coordinates=(0.0, AXLE_DIAMETER, AXLE_LENGTH)), point3=v.findAt(coordinates=(
        0.0, AXLE_DIAMETER, 0.0)))
    p = mdb.models['Model-1'].parts['Connector_beam']
    v = p.vertices
    #Datum plane on bottom inside side face of connector beam slot parallel to arm
    p.DatumPlaneByThreePoints(point1=v.findAt(coordinates=(-AXLE_DIAMETER/4 + OFFSET, AXLE_DIAMETER/4 - OFFSET, 
        AXLE_LENGTH/2 + BASE_WIDTH/2)), point2=v.findAt(coordinates=(-AXLE_DIAMETER/4 + OFFSET, AXLE_DIAMETER/4 - OFFSET, 
        AXLE_LENGTH/2 - BASE_WIDTH/2)), 
        point3=v.findAt(coordinates=(AXLE_DIAMETER/4 + OFFSET, 3*AXLE_DIAMETER/4 - OFFSET, 
        AXLE_LENGTH/2 - BASE_WIDTH/2)))
    #Datum plane on side wall of inside slot on left side if looking at catapult from the rear 
    p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=AXLE_LENGTH/2 + BASE_WIDTH/2)
    connectorBeamDatumIndices['Crossbeam_base_surface_assembly_a2'] = d.keys()[-1]

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
    print('Partitioning side wall1')
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


    print('Partitioning clevis')
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
    #Define Surfaces
    ##################################  
    print('Defining Surfaces')
    ########### Start of Temp Code #################

    ###########creating surfaces on ARM
    # creating surface of base of the arm

    a = mdb.models['Model-1'].rootAssembly
    p = mdb.models['Model-1'].parts['ARM']
    s = p.faces
    side1Faces = s.findAt(((0.01, 0, 0), ), ((0.01, BASE_HEIGHT/2, THICKNESS/2), ), ((0.01, BASE_HEIGHT/2, BASE_WIDTH/2), ), ((0.01, BASE_HEIGHT/2, BASE_WIDTH-THICKNESS/2), 
        ), ((0.01, 0, BASE_WIDTH), ), ((0.01, -BASE_HEIGHT/2, THICKNESS/2), ), ((0.01, -BASE_HEIGHT/2, BASE_WIDTH-THICKNESS/2), ), ((0.01, -BASE_HEIGHT/2, BASE_WIDTH/2), ))
    p.Surface(side1Faces=side1Faces, name='ArmBase')
    #creating surface of the bottom of the spoon on arm
    xOff = SPOON_LENGTH/4
    zOff = BASE_WIDTH/4
    s = p.faces
    side1Faces = s.findAt(((BASE_LENGTH+ARM_LENGTH+SPOON_LENGTH/2+xOff, -TIP_HEIGHT/2, BASE_WIDTH/2+zOff), ), ((BASE_LENGTH+ARM_LENGTH+SPOON_LENGTH/2+xOff, -TIP_HEIGHT/2, BASE_WIDTH/2 -zOff), ), ((BASE_LENGTH+ARM_LENGTH+SPOON_LENGTH/2-xOff, -TIP_HEIGHT/2, BASE_WIDTH/2+zOff), ), ((BASE_LENGTH+ARM_LENGTH+SPOON_LENGTH/2-xOff, -TIP_HEIGHT/2, BASE_WIDTH/2-zOff), ))
    p.Surface(side1Faces=side1Faces, name='SPOON_PULL_SURF')
    
    s = p.faces
    side1Faces = s.findAt(((BASE_LENGTH+ARM_LENGTH+SPOON_LENGTH/2+xOff, TIP_HEIGHT/2, BASE_WIDTH/2+zOff), ), ((BASE_LENGTH+ARM_LENGTH+SPOON_LENGTH/2+xOff, TIP_HEIGHT/2, BASE_WIDTH/2 -zOff), ), ((BASE_LENGTH+ARM_LENGTH+SPOON_LENGTH/2-xOff, TIP_HEIGHT/2, BASE_WIDTH/2+zOff), ), ((BASE_LENGTH+ARM_LENGTH+SPOON_LENGTH/2-xOff, TIP_HEIGHT/2, BASE_WIDTH/2-zOff), ))
    p.Surface(side1Faces=side1Faces, name='SPOON_TOP_SURF')
    #creating surface of the two L members that form spoon 

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

    #Creating Edges on Payload for contact 
    p = mdb.models['Model-1'].parts['Payload']
    e = p.edges
    edges = e.findAt(((0.0, 0.0, PAY_WIDTH/2), ))
    p.Set(edges=edges, name='FRONT_PAYLOAD')
    p = mdb.models['Model-1'].parts['Payload']
    e = p.edges
    edges = e.findAt(((PAY_LENGTH, 0.0, PAY_WIDTH/2), ))
    p.Set(edges=edges, name='BACK_PAYLOAD')
    p = mdb.models['Model-1'].parts['Payload']
    p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=0.0)
    payloadDatumIndices['Pay_Principal_xy'] = p.datums.keys()[-1]
    p = mdb.models['Model-1'].parts['Payload']
    p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=0.0)
    payloadDatumIndices['Pay_FrontPlane_yz'] = p.datums.keys()[-1]
    p = mdb.models['Model-1'].parts['Payload']
    p.DatumPlaneByPrincipalPlane(principalPlane=XZPLANE, offset=0.0)
    payloadDatumIndices['Pay_BottomPlane_xz'] = p.datums.keys()[-1]



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
    a = mdb.models['Model-1'].rootAssembly
    p = mdb.models['Model-1'].parts['CROSSMEMBER']
    a.Instance(name='CROSSMEMBER-1', part=p, dependent=ON)
    p1 = a.instances['CROSSMEMBER-1']
    d1 = a.instances['CROSSMEMBER-1'].datums
    d2 = a.instances['Side_wall_1-1'].datums
    a.FaceToFace(movablePlane=d1[9], fixedPlane=d2[5], flip=OFF, clearance=0.0)
    d1 = a1.instances['CROSSMEMBER-1'].datums
    d2 = a1.instances['Side_wall_2-1'].datums
    a1.CoincidentPoint(movablePoint=d1[10], fixedPoint=d2[4])

    
    #Instance bt ARM and AXLE
    p = mdb.models['Model-1'].parts['ARM']
    a.Instance(name='ARM-1', part=p, dependent=ON)
    a = mdb.models['Model-1'].rootAssembly
    d = p.datums
    d11 = a.instances['ARM-1'].datums
    d12 = a.instances['Connector_beam-1'].datums

    #Constrain arm base to slot base
    a.FaceToFace(movablePlane=d11[armDatumIndices['Arm_base_assembly_plane_yz_a1']], fixedPlane=d12[connectorBeamDatumIndices['Crossbeam_base_surface_assembly_a1']], flip=OFF, clearance=0.0)
    p = mdb.models['Model-1'].parts['Side_wall_1']
    d11 = a.instances['ARM-1'].datums
    d12 = a.instances['Connector_beam-1'].datums
    #Constrain arm side to slot side
    a.FaceToFace(movablePlane=d11[armDatumIndices['Arm_side_assembly_plane_xy_a2']], fixedPlane=d12[connectorBeamDatumIndices['Crossbeam_base_surface_assembly_a2']], flip=OFF, clearance=0.0)
    d11 = a.instances['ARM-1'].datums
    d12 = a.instances['Connector_beam-1'].datums
    #Constrain arm top to slot inside top
    
    d11 = a.instances['ARM-1'].datums
    d12 = a.instances['Connector_beam-1'].datums
    a.FaceToFace(movablePlane=d11[armDatumIndices['Arm_top_assembly_plane_xz_a3']], fixedPlane=d12[connectorBeamDatumIndices['Crossbeam_top_inside_surface_assembly_a3']], flip=ON, clearance=0.0)
    
    
    
    # #################################
    # # Payload Stuff #################
    # #################################
    #Partitioning Payload
    p = mdb.models['Model-1'].parts['Payload']
    c = p.cells
    d = p.datums
    pickedCells = c.findAt(((0, 0, 0), ))
    e, v, d = p.edges, p.vertices, p.datums
    p.PartitionCellByPlanePointNormal(normal=e.findAt(coordinates=(PAY_LENGTH/2, PAY_HEIGHT, 
        PAY_WIDTH)), cells=pickedCells, point=p.InterestingPoint(edge=e.findAt(
        coordinates=(PAY_LENGTH/4, PAY_HEIGHT, 
        PAY_WIDTH)), rule=MIDDLE))
        
    # xOffset = SPOON_FLANGE_HEIGHT/2
    p = mdb.models['Model-1'].parts['Payload']
    p.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE, offset=(PAY_WIDTH/2))
    payloadDatumIndices['Pay_longitudinal_plane_xy'] = d.keys()[-1]
    c = p.cells
    pickedCells = c.findAt(((PAY_LENGTH/4, PAY_HEIGHT/2, PAY_WIDTH/2), ), ((3*PAY_LENGTH/4, PAY_HEIGHT/2, PAY_WIDTH/2), ))
    d = p.datums
    p.PartitionCellByDatumPlane(datumPlane=d[payloadDatumIndices['Pay_longitudinal_plane_xy']], cells=pickedCells)
    # p = mdb.models['Model-1'].parts['Payload']
    # p.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=(PAY_LENGTH-xOffset))

    # c = p.cells
    # pickedCells = c.findAt(((PAY_LENGTH/4, PAY_HEIGHT/2, PAY_WIDTH/2), ), ((3*PAY_LENGTH/4, PAY_HEIGHT/2, PAY_WIDTH/2), ))
    # d2 = p.datums
    #p.PartitionCellByDatumPlane(datumPlane=d2[11], cells=pickedCells)
    
    v = p.vertices
    verts = v.findAt(((PAY_LENGTH/2, 0.0, PAY_WIDTH/2), ))
    p.Set(vertices=verts, name='PAY_PULL_POINT')
    
    s = p.faces
    side1Faces = s.findAt(((PAY_LENGTH/4, 0.0, PAY_WIDTH/4), ), ((PAY_LENGTH/4, 0.0, 3*PAY_WIDTH/4), ), ((3*PAY_LENGTH/4, 0.0, PAY_WIDTH/4), ), ((3*PAY_LENGTH/4, 0.0, 3*PAY_WIDTH/4), ))
    p.Surface(side1Faces=side1Faces, name='PAY_BOTTOM_SURF')
    
    #create bottom node set so that I can use it to enforce contact
    f = p.faces
    faces = f.findAt(((PAY_LENGTH/4, 0.0, PAY_WIDTH/4), ), ((PAY_LENGTH/4, 0.0, 3*PAY_WIDTH/4), ), 
        ((3*PAY_LENGTH/4, 0.0, PAY_WIDTH/4), ), ((3*PAY_LENGTH/4, 0.0, 3*PAY_WIDTH/4), ))
    p.Set(faces=faces, name='BOTTOM_PAYLOAD')
    
    print(payloadDatumIndices)
    print(armDatumIndices)
    ## Adding in payload instance 
    a = mdb.models['Model-1'].rootAssembly
    p = mdb.models['Model-1'].parts['Payload']
    a.Instance(name='Payload-1', part=p, dependent=ON)
    d11 = a.instances['Payload-1'].datums
    d12 = a.instances['ARM-1'].datums
    a.FaceToFace(movablePlane=d11[payloadDatumIndices['Pay_FrontPlane_yz']], fixedPlane=d12[armDatumIndices['spoon_terminus_yz']], flip=OFF, clearance=0.0)
    d11 = a.instances['Payload-1'].datums
    d12 = a.instances['ARM-1'].datums
    a.FaceToFace(movablePlane=d11[payloadDatumIndices['Pay_Principal_xy']], fixedPlane=d12[armDatumIndices['Arm_side_principal_plane_xy']], flip=OFF, clearance=0.0)
    d11 = a.instances['Payload-1'].datums
    d12 = a.instances['ARM-1'].datums
    a.FaceToFace(movablePlane=d11[payloadDatumIndices['Pay_BottomPlane_xz']], fixedPlane=d12[armDatumIndices['horizontal_long_spoon_xz_cut']], flip=OFF, clearance=0.0)
    
    ###############################
    #CLEVIS ROTATION 
    ###############################
    
    #extract locations of pull points
    # 'AXLE_CENTER_POINT' CHANGES RMW 12/3
    spoonCoords = mdb.models['Model-1'].rootAssembly.instances['Payload-1'].sets['PAY_PULL_POINT'].vertices[0].pointOn[0]
    clevisAxleCoords = mdb.models['Model-1'].rootAssembly.instances['Side_wall_1-1'].datums[4].pointOn #NOTE!!! add a dynamic refence to the datum index (like I did for the arm)

    print(spoonCoords)
    print(clevisAxleCoords)

    # spoonCoords = mdb.models['Model-1'].rootAssembly.instances['ARM-1'].sets['SPOON_PULL_POINT'].nodes[0].coordinates



    #Rotating the CLEVIS
    CLEVIS_ANGLE = 0.0 if abs(spoonCoords[0]-clevisAxleCoords[0]) < 0.001 else atan2(spoonCoords[0]-clevisAxleCoords[0],spoonCoords[1]-clevisAxleCoords[1])*180/pi
    print("Updated clevis angle: "+str(CLEVIS_ANGLE))
    a = mdb.models['Model-1'].rootAssembly
    a1.rotate(instanceList=('CROSSMEMBER-1', ), axisPoint=(WALL_LENGTH - WALL_SIDE_X - CLEVIS_RAD, WALL_HEIGHT/2, WALL_THICKNESS/2), 
        axisDirection=(0.0, 0.0, -AXLE_LENGTH), angle=CLEVIS_ANGLE)

    pullPointCoords = mdb.models['Model-1'].rootAssembly.instances['CROSSMEMBER-1'].sets['XBEAM_PULL_POINT'].vertices[0].pointOn[0]
    print(pullPointCoords)
    ropeDistance = sqrt((pullPointCoords[0]-spoonCoords[0])**2 + (pullPointCoords[1]-spoonCoords[1])**2 + (pullPointCoords[2]-spoonCoords[2])**2)
    print(ropeDistance)
    #END CHANGES RMW 12/3

    ###############################
    #WIRE CREATION AND MANIPULATION
    ###############################
    s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
        sheetSize=200.0)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=STANDALONE)
    s.Line(point1=(pullPointCoords[0], pullPointCoords[1]), point2=(spoonCoords[0],spoonCoords[1]))
    p = mdb.models['Model-1'].Part(name='WIRE_BOI', dimensionality=THREE_D, 
        type=DEFORMABLE_BODY)
    p = mdb.models['Model-1'].parts['WIRE_BOI']
    p.BaseWire(sketch=s)
    s.unsetPrimaryObject()
    p = mdb.models['Model-1'].parts['WIRE_BOI']
    session.viewports['Viewport: 1'].setValues(displayedObject=p)
    del mdb.models['Model-1'].sketches['__profile__']

    #create wire material with negative thermal expansion
    CTE_VAL = WIRE_DRAW_DISTANCE*ropeDistance
    
    mdb.models['Model-1'].Material(name='WireMatl')
    mdb.models['Model-1'].materials['WireMatl'].Expansion(table=((-CTE_VAL, ), ))
    mdb.models['Model-1'].materials['WireMatl'].Density(table=((9000.0, ), ))
    mdb.models['Model-1'].materials['WireMatl'].Elastic(table=((10000000000.0, 
        0.3), ))

    mdb.models['Model-1'].HomogeneousSolidSection(name='WireSection', 
        material='WireMatl', thickness=None)
    del mdb.models['Model-1'].sections['WireSection']
    mdb.models['Model-1'].CircularProfile(name='WireCircProfile', r=DRAW_WIRE_RADIUS)
    mdb.models['Model-1'].BeamSection(name='WireSection', 
        integration=DURING_ANALYSIS, poissonRatio=0.0, profile='WireCircProfile', 
        material='WireMatl', temperatureVar=LINEAR, consistentMassMatrix=False)
    p = mdb.models['Model-1'].parts['WIRE_BOI']
    e = p.edges
    # edges = e.findAt(coordinates=((AXLE_LENGTH-WALL_SEP_LENGTH)/2+WALL_SEP_LENGTH/2-XBEAM_FLAT_LENGTH/2, XBEAM_HEIGHT, XBEAM_THICKNESS/4))
    allEdges = e[:]
    region = p.Set(edges=allEdges, name='All_Wire_Edges')
    p = mdb.models['Model-1'].parts['WIRE_BOI']
    p.SectionAssignment(region=region, sectionName='WireSection', offset=0.0, 
        offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)
        
    region=p.sets['All_Wire_Edges']
    p.assignBeamSectionOrientation(region=region, method=N1_COSINES, n1=(0.0, 0.0, 
        -1.0))
    a = mdb.models['Model-1'].rootAssembly
    a.regenerate()

    #Mesh the wire with a very coarse mesh
    p.seedPart(size=1.0, deviationFactor=0.1, minSizeFactor=0.75)
    elemType1 = mesh.ElemType(elemCode=B31, elemLibrary=STANDARD)
    pickedRegions =(allEdges, )
    p.setElementType(regions=pickedRegions, elemTypes=(elemType1, ))
    p.generateMesh()

    #create wire instance to add to assy
    a.Instance(name='WIRE_BOI-1', part=p, dependent=ON)

    #create some geo sets for later
    v1 = a.instances['WIRE_BOI-1'].vertices
    verts1 = v1.findAt(((spoonCoords[0],spoonCoords[1], 0.0), ))
    a.Set(vertices=verts1, name='WireEnd_Arm')
    verts1 = v1.findAt(((pullPointCoords[0],pullPointCoords[1], 0.0), ))
    a.Set(vertices=verts1, name='WireEnd_Clevis')

    #translate instance
    a1.translate(instanceList=('WIRE_BOI-1', ), vector=(0.0, 0.0, spoonCoords[2]))


    ##################################  
    #Define Steps
    ##################################  
    print('Defining the Steps')

    #configure predefined temperature field on the wire
    region = mdb.models['Model-1'].rootAssembly.instances['WIRE_BOI-1'].sets['All_Wire_Edges']
    mdb.models['Model-1'].Temperature(name='WireTemp', createStepName='Initial', 
        region=region, distributionType=UNIFORM, 
        crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, magnitudes=(0.0, ))

    # Defining Loading step
    # mdb.models['Model-1'].ImplicitDynamicsStep(name='Loading', previous='Initial', 
        # maxNumInc=30, application=QUASI_STATIC, initialInc=0.1, nohaf=OFF, 
        # amplitude=RAMP, alpha=DEFAULT, initialConditions=OFF, nlgeom=ON)
    # mdb.models['Model-1'].steps['Loading'].setValues(timePeriod=1.0, maxNumInc=500, 
        # initialInc=0.1, minInc=5e-05)
    mdb.models['Model-1'].StaticStep(name='Loading', previous='Initial', 
        maxNumInc=30, initialInc=0.1, nlgeom=ON)
    mdb.models['Model-1'].steps['Loading'].setValues(timePeriod=1.0, maxNumInc=500, minInc=5e-05)
        
    # Defining Launch step
    mdb.models['Model-1'].ImplicitDynamicsStep(name='Launch', previous='Loading', 
        timePeriod=0.01, application=TRANSIENT_FIDELITY, initialInc=0.01, 
        minInc=2e-07, nohaf=OFF, initialConditions=ON)
    mdb.models['Model-1'].steps['Launch'].setValues(maxNumInc=500)
    
    # Defining Buckling step
    mdb.models['Model-1'].BuckleStep(name='BuckleCheck', previous='Loading', numEigen=1, 
        eigensolver=LANCZOS, minEigen=None, blockSize=DEFAULT, maxBlocks=DEFAULT)
    
    # Defining FollowThru step
    mdb.models['Model-1'].ImplicitDynamicsStep(name='FollowThru', previous='Launch', 
        timePeriod=0.25, application=TRANSIENT_FIDELITY, initialInc=0.015, 
        minInc=2e-05, nohaf=OFF, initialConditions=ON)
    mdb.models['Model-1'].steps['FollowThru'].setValues(maxNumInc=300) ####JOB BROKE WITH ONLY 100 INCS LAST TIME   
       
    #update save increments for ODB output
    mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValuesInStep(
        stepName='Loading', timeInterval=0.1)
    mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValuesInStep(
        stepName='Launch', timeInterval=0.005)
    mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValuesInStep(
        stepName='FollowThru', timeInterval=0.025)

    #  approximate timesteps save to ODB
    mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(
        timeMarks=OFF)

    mdb.models['Model-1'].SmoothStepAmplitude(name='Amp-1', timeSpan=STEP, data=((
        0.0, 0.0), (1.0, 1.0)))
    mdb.models['Model-1'].predefinedFields['WireTemp'].setValuesInStep(
        stepName='Loading', magnitudes=(1.0, ), amplitude='Amp-1')

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
    region1=a.instances['CROSSMEMBER-1'].surfaces['RightSide']
    region2=a.instances['Side_wall_1-1'].surfaces['CircularNotch']
    mdb.models['Model-1'].Tie(name='Cross__Sidewall1', master=region1, 
        slave=region2, positionToleranceMethod=COMPUTED, adjust=ON, 
        tieRotations=ON, thickness=ON)
    region1=a.instances['CROSSMEMBER-1'].surfaces['LeftSide']
    region2=a.instances['Side_wall_2-1'].surfaces['CircularNotch']
    mdb.models['Model-1'].Tie(name='Cross__Sidewall2', master=region1, 
        slave=region2, positionToleranceMethod=COMPUTED, adjust=ON, 
        tieRotations=ON, thickness=ON)
        
    #Tie Wire to the clevis
    region1=a.instances['CROSSMEMBER-1'].surfaces['CLEVIS_TOP']
    region2=a.sets['WireEnd_Clevis']
    mdb.models['Model-1'].Tie(name='Lower_Wire_Tie', master=region1, slave=region2, 
        positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)
        
    # #Tie Wire to the Spoon Back
    # region1=a.instances['ARM-1'].surfaces['SPOON_PULL_SURF']
    # region2=a.sets['WireEnd_Arm']
    # mdb.models['Model-1'].Tie(name='Upper_Wire_Tie', master=region1, slave=region2, 
        # positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)
    
    #Tie Wire to the Payload Bottom
    region1=a.instances['Payload-1'].surfaces['PAY_BOTTOM_SURF']
    region2=a.sets['WireEnd_Arm']
    mdb.models['Model-1'].Tie(name='Upper_Wire_Tie', master=region1, slave=region2, 
        positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)
        
    #Wire deletion step
    a = mdb.models['Model-1'].rootAssembly
    edges1 = region = mdb.models['Model-1'].rootAssembly.instances['WIRE_BOI-1'].sets['All_Wire_Edges']
    mdb.models['Model-1'].ModelChange(name='WireDeletion', createStepName='Launch', 
        region=region, activeInStep=False, includeStrain=False)
    #: The interaction "Int-2" has been created.
    # mdb.models['Model-1'].interactions['Int-1'].reset('Launch')

    #ADD CONTACT
    #define contact property
    mdb.models['Model-1'].ContactProperty('PayloadRoughContactProp')
    mdb.models['Model-1'].interactionProperties['PayloadRoughContactProp'].TangentialBehavior(
        formulation=ROUGH)
    mdb.models['Model-1'].interactionProperties['PayloadRoughContactProp'].NormalBehavior(
        pressureOverclosure=LINEAR, contactStiffness=10000000000.0, 
        constraintEnforcementMethod=DEFAULT)
    
    #Create non-sliding contact b/n payload and arm
    region1=a.instances['ARM-1'].surfaces['SPOON_TOP_SURF']
    region2=a.instances['Payload-1'].surfaces['PAY_BOTTOM_SURF']
    regionDef=mdb.models['Model-1'].rootAssembly.allInstances['Payload-1'].sets['BOTTOM_PAYLOAD']
    mdb.models['Model-1'].SurfaceToSurfaceContactStd(name='ROUGH_CONTACT', 
        createStepName='Initial', master=region1, slave=region2, sliding=SMALL, 
        thickness=ON, interactionProperty='PayloadRoughContactProp', 
        adjustMethod=SET, initialClearance=OMIT, datumAxis=None, 
        clearanceRegion=None, tied=OFF, adjustSet=regionDef)

    ##################################  
    #Create Loads
    ##################################  
    print('Defining Loads')
    #creating gravity
    mdb.models['Model-1'].Gravity(name='Load-1', createStepName='Loading', 
        comp2=-1.62, distributionType=UNIFORM, field='')
        
    #BuckleCheck Dummy Load
    region = a.instances['Payload-1'].sets['PAY_PULL_POINT']
    mdb.models['Model-1'].ConcentratedForce(name='BuckleDummy', 
        createStepName='BuckleCheck', region=region, cf1=1.0, cf2=-1.0,
        distributionType=UNIFORM, field='', localCsys=None)   
        
        
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

    #prevent wire from rotating about its axis due to roundoff error.
    region = mdb.models['Model-1'].rootAssembly.instances['WIRE_BOI-1'].sets['All_Wire_Edges']
    mdb.models['Model-1'].DisplacementBC(name='WireNoVertRotation', createStepName='Initial', 
        region=region, u1=UNSET, u2=UNSET, u3=UNSET, ur1=SET, ur2=SET, ur3=UNSET, 
        amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)

    #prevent payload from sliding in Z direction
    # region = a.instances['Payload-1'].sets['Set-1']
    # mdb.models['Model-1'].DisplacementBC(name='PayloadNoZSlide', 
        # createStepName='Initial', region=region, u1=UNSET, u2=UNSET, u3=SET, 
        # ur1=UNSET, ur2=UNSET, ur3=UNSET, amplitude=UNSET, distributionType=UNIFORM, 
        # fieldName='', localCsys=None)
    ##################################      
    #Mesh Parts
    ##################################  
    print('Meshing the Part')

    ##Meshing Arm
    p = mdb.models['Model-1'].parts['ARM']
    p.seedPart(size=0.1, deviationFactor=0.1, minSizeFactor=0.1)
    p.generateMesh()

    ##Meshing Crossmember
    p = mdb.models['Model-1'].parts['CROSSMEMBER']
    c = p.cells
    pickedRegions = c.findAt(((0,0,0), ), ((0,0,AXLE_DIAMETER/2), ), ((AXLE_LENGTH/2-0.001,0,0), ), ((AXLE_LENGTH/2-0.001, 0, 
        AXLE_DIAMETER/2), ), ((AXLE_LENGTH/2+0.001,0,0), ), ((AXLE_LENGTH/2+0.001, 0, 
        AXLE_DIAMETER/2), ), ((AXLE_LENGTH,0,0), ), ((AXLE_LENGTH,0,AXLE_DIAMETER/2), ))
    p.setMeshControls(regions=pickedRegions, elemShape=TET, technique=FREE)
    elemType1 = mesh.ElemType(elemCode=C3D20R)
    elemType2 = mesh.ElemType(elemCode=C3D15)
    elemType3 = mesh.ElemType(elemCode=C3D10)
    pickedRegions =(cells, )
    p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, 
        elemType3))
    p.seedPart(size=0.025, deviationFactor=0.1, minSizeFactor=0.1)
    p.generateMesh()

    ##Meshing Connector beam
    p = mdb.models['Model-1'].parts['Connector_beam']
    p.seedPart(size=0.036, deviationFactor=0.1, minSizeFactor=0.1)
    p.generateMesh()

    ##Meshing Payload
    p = mdb.models['Model-1'].parts['Payload']
    p.seedPart(size=0.019, deviationFactor=0.1, minSizeFactor=0.1)
    p.generateMesh()

    ##Meshing Side wall 1
    p = mdb.models['Model-1'].parts['Side_wall_1']
    p.seedPart(size=0.04, deviationFactor=0.1, minSizeFactor=0.1)
    p.generateMesh()

    ##Meshing Side wall 2
    p = mdb.models['Model-1'].parts['Side_wall_2']
    p.seedPart(size=0.04, deviationFactor=0.1, minSizeFactor=0.1)
    p.generateMesh()

    ##Meshing wire
    p = mdb.models['Model-1'].parts['WIRE_BOI']
    session.viewports['Viewport: 1'].setValues(displayedObject=p)
    p.generateMesh()
    
    ####################################
    ### CREATE ELEMENT SET (EXCLUDE WIRE AND PAYLOAD
    ####################################
    a = mdb.models['Model-1'].rootAssembly
    a.regenerate()
    elements1 = a.instances['Side_wall_1-1'].elements[:]
    elements2 = a.instances['Connector_beam-1'].elements[:]
    elements3 = a.instances['Side_wall_2-1'].elements[:]
    elements4 = a.instances['CROSSMEMBER-1'].elements[:]
    elements5 = a.instances['ARM-1'].elements[:]
    
    a.Set(elements=elements1+elements2+elements3+elements4+elements5, 
        name='ALL_PART')

    # ####################################
    # ## Creation/Execution of the Job ###
    # ####################################
    print 'Creating/Running Job'
    ModelName='Model-1'
    JobName = 'Job-1'
    mdb.Job(name=JobName, model=ModelName, description='', type=ANALYSIS, 
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
        scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=4, 
        numDomains=4, numGPUs=0)

    job=mdb.jobs[JobName]

    #delete lock file, which for some reason tends to hang around, if it exists
    if os.access('%s.lck'%JobName,os.F_OK):
        os.remove('%s.lck'%JobName)
           
    job.submit()
    job.waitForCompletion()
    print 'Completed job'

    ##END LOOP (i.e., end indentation)

    ##################################      
    #Output Variables
    ################################## 

    # Mass
    prop = mdb.models[ModelName].rootAssembly.getMassProperties()
    mass = prop['mass']

    # # Payload exit velocity & angle
    # velocity_max, veloAngle = Post_P_Script_Velo.getResults(JobName)

    # Max Mises stress in structure
    outputList = Post_P_Script_Velo.getResults(JobName) 
    
    return outputList
