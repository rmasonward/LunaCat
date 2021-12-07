"""
...
"""

from abaqus import *
from abaqusConstants import *
import visualization
from viewerModules import *

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def getResults(ModelName):

    """
    This ODB reading script does the following:
    -Retrieves the displacement at TIPNODE
    -Scans for max. Mises stress in a part (if set exists)
    """
    
    # Open the output database.
    print 'Made it in'
    odbName = ModelName+'.odb'
    odb = visualization.openOdb(odbName)    
   
    print 'Scanning the PART for maximum VM STRESS'
    elsetName='ALL_PART'
    elset = elemset = None
    region = "over the entire model"
    assembly = odb.rootAssembly

    #Check to see if the element set exists
    #in the assembly

    if elsetName:
       try:
           elemset = assembly.elementSets[elsetName]
           region = " in the element set : " + elsetName;
       except KeyError:
           print 'An assembly level elset named %s does' \
                  'not exist in the output database %s' \
                  % (elsetName, odbName)
           odb.close()
           exit(0)
           
    """ Initialize maximum values """
    maxMises = -0.1
    maxVMElem = 0
    maxStep = "_None_"
    maxFrame = -1
    Stress = 'S'
    isStressPresent = 0
    for step in odb.steps.values():
       print 'Processing Step:', step.name
       for frame in step.frames:
           allFields = frame.fieldOutputs
           if (allFields.has_key(Stress)):
               isStressPresent = 1
               stressSet = allFields[Stress]
               if elemset:
                   stressSet = stressSet.getSubset(
                       region=elemset)      
               for stressValue in stressSet.values:                
                   if (stressValue.mises > maxMises):
                       maxMises = stressValue.mises
                       maxVMElem = stressValue.elementLabel
                       maxStep = step.name
                       maxFrame = frame.incrementNumber
    if(isStressPresent):
       print 'Maximum von Mises stress %s is %f in element %d'%(
           region, maxMises, maxVMElem)
       print 'Location: frame # %d  step:  %s '%(maxFrame,maxStep)
    else:
       print 'Stress output is not available in' \
             'the output database : %s\n' %(odb.name)     

    odb.close()

    return maxMises


   
