"""
...
"""

from abaqus import *
from abaqusConstants import *
import visualization
from viewerModules import *
import math

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def getResults(ModelName):

    """
    This ODB reading script does the following:
    -Scans to see when the there is no more contact between the payload and arm
    -Pulls the velocity of the payload at this point
    """
    
    # Open the output database.
    print 'Made it in'
    odbName = ModelName+'.odb'
    print odbName
    odb = visualization.openOdb(odbName)
    # lastFrame = odb.steps['Step-1'].frames[-1]
    launchFrames = odb.steps['FollowThru'].frames
    print 'odb open'
 

    # Selecting the node(s) to be queried
    pTip = odb.rootAssembly.instances['ARM-1'].nodeSets['SPOON_PULL_POINT']
    payPoint = odb.rootAssembly.instances['Payload-1'].nodeSets['PAY_PULL_POINT']
    
    veloMagMax = 0
    veloAngle = 0
    
    print("Found pTip")
    
        #looking for when there is no contact and taking velocity at that point
    for frame in launchFrames:
        framecontact = frame.fieldOutputs['CPRESS'].getSubset(region=payPoint).values[0].data
        #want V1, V2
        if framecontact < 0.0001:
            frameVelo = frame.fieldOutputs['V'].getSubset(region=pTip).values[0].data

            veloX = frameVelo[0]
            veloY = frameVelo[1]
            
            frameVeloMag = math.sqrt(veloX**2 + veloY**2)
            if frameVeloMag > veloMagMax:
                veloMagMax = frameVeloMag
                veloAngle = math.atan2(veloX,veloY)
    
    print("Found the maximum velocity and angle")
    
    
    # # Retrieve Y-displacements at the splines/connectors
    # print 'Retrieving ALL final displacements at ALL points'
    # dispField = lastFrame.fieldOutputs['U']

    # print 'Retrieving ALL displacements at TIPNODE'
    # dFieldpTip = dispField.getSubset(region=pTip)
    
    # print 'Retrieving only U2 at TIPNODE'
    # #Note, U1=data[0], U2=data[1], U3=data[2]
    # disppTip = dFieldpTip.values[0].data[1]
    
   
	## The following is left for use in later probems/projects
##    print 'Scanning the PART for maximum VM STRESS'
##    elsetName='ALL_PART'
##    elset = elemset = None
##    region = "over the entire model"
##    assembly = odb.rootAssembly
##
##    #Check to see if the element set exists
##    #in the assembly
##    
##    if elsetName:
##        try:
##            elemset = assembly.elementSets[elsetName]
##            region = " in the element set : " + elsetName;
##        except KeyError:
##            print 'An assembly level elset named %s does' \
##                   'not exist in the output database %s' \
##                   % (elsetName, odbName)
##            odb.close()
##            exit(0)
##            
##    """ Initialize maximum values """
##    maxMises = -0.1
##    maxVMElem = 0
##    maxStep = "_None_"
##    maxFrame = -1
##    Stress = 'S'
##    isStressPresent = 0
##    for step in odb.steps.values():
##        print 'Processing Step:', step.name
##        for frame in step.frames:
##            allFields = frame.fieldOutputs
##            if (allFields.has_key(Stress)):
##                isStressPresent = 1
##                stressSet = allFields[Stress]
##                if elemset:
##                    stressSet = stressSet.getSubset(
##                        region=elemset)      
##                for stressValue in stressSet.values:                
##                    if (stressValue.mises > maxMises):
##                        maxMises = stressValue.mises
##                        maxVMElem = stressValue.elementLabel
##                        maxStep = step.name
##                        maxFrame = frame.incrementNumber
##    if(isStressPresent):
##        print 'Maximum von Mises stress %s is %f in element %d'%(
##            region, maxMises, maxVMElem)
##        print 'Location: frame # %d  step:  %s '%(maxFrame,maxStep)
##    else:
##        print 'Stress output is not available in' \
##              'the output database : %s\n' %(odb.name)     

    odb.close()
    
    return veloMagMax,veloAngle

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def createXYPlot(vpOrigin, vpName, plotName, data):

	##NOTE: I have never used this but it might be interesting in some problems
    
    """
    Display curves of theoretical and computed results in
    a new viewport.
    """
    print 'In plotter'
    from visualization import  USER_DEFINED
    
    vp = session.Viewport(name=vpName, origin=vpOrigin, 
        width=150, height=100)
    print 'Viewport created'
    xyPlot = session.XYPlot(plotName)
    chart = xyPlot.charts.values()[0]
    curveList = []
    for elemName, xyValues in sorted(data.items()):
        xyData = session.XYData(elemName, xyValues)
        curve = session.Curve(xyData)
        curveList.append(curve)

    print 'Curves created'
    chart.setValues(curvesToPlot=curveList)
    chart.axes1[0].axisData.setValues(useSystemTitle=False,title='Arc Height')
    chart.axes2[0].axisData.setValues(useSystemTitle=False,title=plotName)
    vp.setValues(displayedObject=xyPlot)
    print 'Plot displayed'
    return
