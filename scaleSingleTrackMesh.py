''' This codes parses abaqus .inp (reference) file, translates the given nodal coordinates first,
    then scales nodal coordinates of different sets to adjust the deposit length, width, depth and 
    base plate length, width, depth as per user supplied scaling factors and outputs the new nodal 
    coordinates to a file. 
    
    This code avoids repeated creation of single track FE mesh each time the deposit or base plate 
    dimensions change.
    
    NOTE: This code requires that the original nodes be numbered starting from 1.
    If this is not the case, then the code 'renumberNodesAndConnectivity' must 
    be executed prior to using this code
'''
import sys
import os
import statistics as stats
import random as rand
import time
import numpy as np
import math
import inspect
from parseAbqInpFileForNodalCoords import parseAbqInpFileForNodalCoords
from parseAbqInpFileForNodeSets import parseAbqInpFileForNodeSets


def main():
    
    nArguments = len(sys.argv)
    if nArguments != 16:
        text = """15 command line arguments are expected: \
                \n\t(1)  ABQ inp file name including path \
                \n\t ****** -------- ******* \
                \n\t(2)  BasePlateLength-left-coarseRegion (original = 3 units) \
                \n\t(3)  BasePlateLength-left-fineRegion (original = 2 units) \
                \n\t(4)  DepositLength (original = 10 units) \
                \n\t(5)  BasePlateLength-right-fineRegion (original = 2 units) \
                \n\t(6)  BasePlateLength-right-coarseRegion (original = 3 units) \
                \n\t ****** -------- ******* \
                \n\t(7)  BasePlateWidth-left-coarseRegion (original = 1 unit) \
                \n\t(8)  BasePlateWidth-left-fineRegion (original = 1 unit) \
                \n\t(9)  DepositWidth (original = 1 unit) \
                \n\t(10) BasePlateWidth-right-fineRegion (original = 1 unit) \
                \n\t(11) BasePlateWidth-right-coarseRegion (original = 1 unit) \
                \n\t ****** -------- ******* \
                \n\t(12) BasePlateDepth-bottom-coarseRegion (original = 1 unit) \
                \n\t(13) BasePlateDepth-top-fineRegion (original = 0.75 units) \
                \n\t(14) DepositDepth (original = 0.5 unit) \
                \n\t ****** -------- ******* \
                \n\t(15) Output directory path"""
        raise RuntimeError(text)

    nameFileInpAbaqus = sys.argv[1]
    scalingFactorBasePlateLengthLeftCoarseRegion = float(sys.argv[2])/3 # Divide target dimensions with reference dimensions to obtain scaling factors
    scalingFactorBasePlateLengthLeftFineRegion = float(sys.argv[3])/2
    scalingFactorDepositLength = float(sys.argv[4])/10
    scalingFactorBasePlateLengthRightFineRegion = float(sys.argv[5])/2
    scalingFactorBasePlateLengthRightCoarseRegion = float(sys.argv[6])/3
    scalingFactorBasePlateWidthLeftCoarseRegion = float(sys.argv[7])/1
    scalingFactorBasePlateWidthLeftFineRegion = float(sys.argv[8])/1
    scalingFactorDepositWidth = float(sys.argv[9])/1
    scalingFactorBasePlateWidthRightFineRegion = float(sys.argv[10])/1
    scalingFactorBasePlateWidthRightCoarseRegion = float(sys.argv[11])/1
    scalingFactorBasePlateDepthBottomCoarseRegion = float(sys.argv[12])/1
    scalingFactorBasePlateDepthTopFineRegion = float(sys.argv[13])/0.75
    scalingFactorDepositDepth = float(sys.argv[14])/0.5
    outputPath = sys.argv[15]
        
    pathDir = outputPath + "/"    
    # -------------------------------------------------------------------------------------------------------------------------------
    # -------------------------------------- Parse the ABQ inp file and obtain nodal coord information
    # -------------------------------------------------------------------------------------------------------------------------------
    listNodes, listCoordinatesX, listCoordinatesY, listCoordinatesZ = parseAbqInpFileForNodalCoords(nameFileInpAbaqus)    
    listNodeSetsNames, listNodesInNodeSets = parseAbqInpFileForNodeSets(nameFileInpAbaqus)
    
    coordXmin = min(listCoordinatesX)
    coordYmin = min(listCoordinatesY)
    coordZmin = min(listCoordinatesZ)
    # Translate the nodal coordinates so that all the coordinates are +ve
    listCoordinatesX = [listCoordinatesX[ii] - coordXmin for ii in range(0,len(listCoordinatesX))]
    listCoordinatesY = [listCoordinatesY[ii] - coordYmin for ii in range(0,len(listCoordinatesY))]
    listCoordinatesZ = [listCoordinatesZ[ii] - coordZmin for ii in range(0,len(listCoordinatesZ))]
    
    #
    listMinCoordinatesForUnorderedSetsAlongLength, listMaxCoordinatesForUnorderedSetsAlongLength = [[] for ii in range(0,2)]
    listMinCoordinatesForUnorderedSetsAlongWidth, listMaxCoordinatesForUnorderedSetsAlongWidth = [[] for ii in range(0,2)]
    listMinCoordinatesForUnorderedSetsAlongDepth, listMaxCoordinatesForUnorderedSetsAlongDepth = [[] for ii in range(0,2)]
    for iNodeSet in range(0,len(listNodeSetsNames)):
        listCoordinatesX_ = [listCoordinatesX[idNode-1] for idNode in listNodesInNodeSets[iNodeSet]]
        listCoordinatesY_ = [listCoordinatesY[idNode-1] for idNode in listNodesInNodeSets[iNodeSet]]
        listCoordinatesZ_ = [listCoordinatesZ[idNode-1] for idNode in listNodesInNodeSets[iNodeSet]]
        
        listMinCoordinatesForUnorderedSetsAlongLength.append(min(listCoordinatesY_))
        listMaxCoordinatesForUnorderedSetsAlongLength.append(max(listCoordinatesY_))

        listMinCoordinatesForUnorderedSetsAlongWidth.append(min(listCoordinatesX_))
        listMaxCoordinatesForUnorderedSetsAlongWidth.append(max(listCoordinatesX_))
       
        listMinCoordinatesForUnorderedSetsAlongDepth.append(min(listCoordinatesZ_))
        listMaxCoordinatesForUnorderedSetsAlongDepth.append(max(listCoordinatesZ_))
    
    # Determine the ordering of node sets in all directions based on increasing coordinate directions
    # Find the ordering of node sets in each direction for applying scaling
    for iNodeSet in range(0,len(listNodeSetsNames)):
        if listNodeSetsNames[iNodeSet].upper() == "setControllingBasePlateLengthLeftMost".upper():
           minYalongLength_left = listMinCoordinatesForUnorderedSetsAlongLength[iNodeSet]
        if listNodeSetsNames[iNodeSet].upper() == "setControllingBasePlateLengthRightMost".upper():
           minYalongLength_right = listMinCoordinatesForUnorderedSetsAlongLength[iNodeSet]

        if listNodeSetsNames[iNodeSet].upper() == "setControllingBasePlateWidthLeftMost".upper():
           minXalongWidth_left = listMinCoordinatesForUnorderedSetsAlongWidth[iNodeSet]
        if listNodeSetsNames[iNodeSet].upper() == "setControllingBasePlateWidthRightMost".upper():
           minXalongWidth_right = listMinCoordinatesForUnorderedSetsAlongWidth[iNodeSet]
    
        if listNodeSetsNames[iNodeSet].upper() == "setControllingBasePlateDepthBottom".upper():
           minZalongDepth_bottom = listMinCoordinatesForUnorderedSetsAlongDepth[iNodeSet]
        if listNodeSetsNames[iNodeSet].upper() == "setControllingDepositDepth".upper():
           listCoordinatesZ_ = [listCoordinatesZ[idNode-1] for idNode in listNodesInNodeSets[iNodeSet]]
           minZalongDepth_top = listMinCoordinatesForUnorderedSetsAlongDepth[iNodeSet]
    
    if minYalongLength_left < minYalongLength_right:
        listNodeSetsForOrderingAlongLength = ['setControllingBasePlateLengthLeftMost',
                                              'setControllingBasePlateLengthLeft',
                                              'setControllingDepositLength',
                                              'setControllingBasePlateLengthRight',
                                              'setControllingBasePlateLengthRightMost']
                                              
        listScalingFactorsAlongLength = [scalingFactorBasePlateLengthLeftCoarseRegion,
                                         scalingFactorBasePlateLengthLeftFineRegion,
                                         scalingFactorDepositLength,
                                         scalingFactorBasePlateLengthRightFineRegion,
                                         scalingFactorBasePlateLengthRightCoarseRegion]
        
    else:
        listNodeSetsForOrderingAlongLength = ['setControllingBasePlateLengthRightMost',
                                              'setControllingBasePlateLengthRight',
                                              'setControllingDepositLength',
                                              'setControllingBasePlateLengthLeft',
                                              'setControllingBasePlateLengthLeftMost']
                                              
        listScalingFactorsAlongLength = [scalingFactorBasePlateLengthRightCoarseRegion,
                                         scalingFactorBasePlateLengthRightFineRegion,
                                         scalingFactorDepositLength,
                                         scalingFactorBasePlateLengthLeftFineRegion,
                                         scalingFactorBasePlateLengthLeftCoarseRegion]
                                         
    listMinCoordinatesForSetOrderingAlongLength = []
    listMaxCoordinatesForSetOrderingAlongLength = []
    for iNodeSetOrdered in range(0,len(listNodeSetsForOrderingAlongLength)):
        for iNodeSet in range(0,len(listNodeSetsNames)):
            if listNodeSetsNames[iNodeSet].upper() == listNodeSetsForOrderingAlongLength[iNodeSetOrdered].upper():
                listMinCoordinatesForSetOrderingAlongLength.append(listMinCoordinatesForUnorderedSetsAlongLength[iNodeSet])
                listMaxCoordinatesForSetOrderingAlongLength.append(listMaxCoordinatesForUnorderedSetsAlongLength[iNodeSet])
                                         
    if minXalongWidth_left < minXalongWidth_right:
        listNodeSetsForOrderingAlongWidth = ['setControllingBasePlateWidthLeftMost',
                                             'setControllingBasePlateWidthLeft',
                                             'setControllingDepositWidth',
                                             'setControllingBasePlateWidthRight',
                                             'setControllingBasePlateWidthRightMost']
                                             
        listScalingFactorsAlongWidth = [scalingFactorBasePlateWidthLeftCoarseRegion,
                                        scalingFactorBasePlateWidthLeftFineRegion,
                                        scalingFactorDepositWidth,
                                        scalingFactorBasePlateWidthRightFineRegion,
                                        scalingFactorBasePlateWidthRightCoarseRegion]
                                        
    else:
        listNodeSetsForOrderingAlongWidth = ['setControllingBasePlateWidthRightMost',
                                             'setControllingBasePlateWidthRight',
                                             'setControllingDepositWidth',
                                             'setControllingBasePlateWidthLeft',
                                             'setControllingBasePlateWidthLeftMost']
                                             
        listScalingFactorsAlongWidth = [scalingFactorBasePlateWidthRightCoarseRegion,
                                        scalingFactorBasePlateWidthRightFineRegion,
                                        scalingFactorDepositWidth,
                                        scalingFactorBasePlateWidthLeftFineRegion,
                                        scalingFactorBasePlateWidthLeftCoarseRegion]

    listMinCoordinatesForSetOrderingAlongWidth = []
    listMaxCoordinatesForSetOrderingAlongWidth = []
    for iNodeSetOrdered in range(0,len(listNodeSetsForOrderingAlongWidth)):
        for iNodeSet in range(0,len(listNodeSetsNames)):
            if listNodeSetsNames[iNodeSet].upper() == listNodeSetsForOrderingAlongWidth[iNodeSetOrdered].upper():
                listMinCoordinatesForSetOrderingAlongWidth.append(listMinCoordinatesForUnorderedSetsAlongWidth[iNodeSet])
                listMaxCoordinatesForSetOrderingAlongWidth.append(listMaxCoordinatesForUnorderedSetsAlongWidth[iNodeSet])
                
    if minZalongDepth_bottom < minZalongDepth_top:
        listNodeSetsForOrderingAlongDepth = ['setControllingBasePlateDepthBottom',
                                             'setControllingBasePlateDepthTop',
                                             'setControllingDepositDepth']
                                             
        listScalingFactorsAlongDepth = [scalingFactorBasePlateDepthBottomCoarseRegion,
                                        scalingFactorBasePlateDepthTopFineRegion,
                                        scalingFactorDepositDepth]
                                        
    else:
        listNodeSetsForOrderingAlongDepth = ['setControllingDepositDepth',
                                             'setControllingBasePlateDepthTop',
                                             'setControllingBasePlateDepthBottom']
                                             
        listScalingFactorsAlongDepth = [scalingFactorDepositDepth,
                                        scalingFactorBasePlateDepthTopFineRegion,
                                        scalingFactorBasePlateDepthBottomCoarseRegion]
                                        
    listMinCoordinatesForSetOrderingAlongDepth = []
    listMaxCoordinatesForSetOrderingAlongDepth = []
    for iNodeSetOrdered in range(0,len(listNodeSetsForOrderingAlongDepth)):
        for iNodeSet in range(0,len(listNodeSetsNames)):
            if listNodeSetsNames[iNodeSet].upper() == listNodeSetsForOrderingAlongDepth[iNodeSetOrdered].upper():
                listMinCoordinatesForSetOrderingAlongDepth.append(listMinCoordinatesForUnorderedSetsAlongDepth[iNodeSet])
                listMaxCoordinatesForSetOrderingAlongDepth.append(listMaxCoordinatesForUnorderedSetsAlongDepth[iNodeSet])                                        
    
    # Scale length dimensions
    translationTotal = 0
    listCoordinatesYscaled = [1E20]*len(listCoordinatesY)
    for iNodeSetOrdering in range(0,len(listNodeSetsForOrderingAlongLength)):
        nameNodeSetOrdering = listNodeSetsForOrderingAlongLength[iNodeSetOrdering]
        minForThisRegion = listMinCoordinatesForSetOrderingAlongLength[iNodeSetOrdering]
        maxForThisRegion = listMaxCoordinatesForSetOrderingAlongLength[iNodeSetOrdering]
        dLength = maxForThisRegion - minForThisRegion     
        for iNodeSet in range(0,len(listNodeSetsNames)):
            if listNodeSetsNames[iNodeSet].upper() == nameNodeSetOrdering.upper():
                for jNodeID in listNodesInNodeSets[iNodeSet]:
                    listCoordinatesYscaled[jNodeID-1] = translationTotal +  minForThisRegion + \
                                    (listCoordinatesY[jNodeID-1] - minForThisRegion)*listScalingFactorsAlongLength[iNodeSetOrdering]
        translationTotal = translationTotal + dLength*(listScalingFactorsAlongLength[iNodeSetOrdering] - 1)
        
    # Scale width dimensions
    translationTotal = 0
    listCoordinatesXscaled = [1E20]*len(listCoordinatesX)
    for iNodeSetOrdering in range(0,len(listNodeSetsForOrderingAlongWidth)):
        nameNodeSetOrdering = listNodeSetsForOrderingAlongWidth[iNodeSetOrdering]
        minForThisRegion = listMinCoordinatesForSetOrderingAlongWidth[iNodeSetOrdering]
        maxForThisRegion = listMaxCoordinatesForSetOrderingAlongWidth[iNodeSetOrdering]
        dLength = maxForThisRegion - minForThisRegion     
        for iNodeSet in range(0,len(listNodeSetsNames)):
            if listNodeSetsNames[iNodeSet].upper() == nameNodeSetOrdering.upper():
                for jNodeID in listNodesInNodeSets[iNodeSet]:
                    listCoordinatesXscaled[jNodeID-1] = translationTotal +  minForThisRegion + \
                                    (listCoordinatesX[jNodeID-1] - minForThisRegion)*listScalingFactorsAlongWidth[iNodeSetOrdering]
        translationTotal = translationTotal + dLength*(listScalingFactorsAlongWidth[iNodeSetOrdering] - 1)
    
    # Scale depth dimensions
    translationTotal = 0
    listCoordinatesZscaled = [1E20]*len(listCoordinatesZ)
    for iNodeSetOrdering in range(0,len(listNodeSetsForOrderingAlongDepth)):
        nameNodeSetOrdering = listNodeSetsForOrderingAlongDepth[iNodeSetOrdering]
        minForThisRegion = listMinCoordinatesForSetOrderingAlongDepth[iNodeSetOrdering]
        maxForThisRegion = listMaxCoordinatesForSetOrderingAlongDepth[iNodeSetOrdering]
        dLength = maxForThisRegion - minForThisRegion     
        for iNodeSet in range(0,len(listNodeSetsNames)):
            if listNodeSetsNames[iNodeSet].upper() == nameNodeSetOrdering.upper():
                for jNodeID in listNodesInNodeSets[iNodeSet]:
                    listCoordinatesZscaled[jNodeID-1] = translationTotal +  minForThisRegion + \
                                    (listCoordinatesZ[jNodeID-1] - minForThisRegion)*listScalingFactorsAlongDepth[iNodeSetOrdering]
        translationTotal = translationTotal + dLength*(listScalingFactorsAlongDepth[iNodeSetOrdering] - 1)    
        
    # write the scaled nodal information to a file
    out_path = pathDir + 'nodalCoordsScaled.dat'
    with open(out_path, 'w') as file_out:
        for iNode in range(0,len(listNodes)):
            file_out.write("{0:8d}{1:2s}{2:25.10f}{3:2s}{4:25.10f}{5:2s}{6:25.10f}\n".format(listNodes[iNode], ",",\
                                                                           listCoordinatesXscaled[iNode], ",",\
                                                                           listCoordinatesYscaled[iNode], ",",\
                                                                           listCoordinatesZscaled[iNode]))
    file_out.close()
        
    
    
if __name__ == "__main__":
    main()