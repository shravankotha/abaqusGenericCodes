''' This codes parses abaqus .inp file and outputs the following
    (1) nodes & nodal coordinates - original
    (2) nodes & nodal coordinates - renumbered to start with number 1
    It also takes element connectivity as input and outputs
    (1) element connectivity - renumbered to account for the above node renumbering 
'''
import sys
import os
import statistics as stats
import random as rand
import time
import numpy as np
import math
from mpi4py import MPI
import inspect
from parseAbqInpFileForNodalCoords import parseAbqInpFileForNodalCoords
from performLinearPartitioning import performLinearPartitioning
from readConnectivityFromAFile import readConnectivity

def main():
    comm = MPI.COMM_WORLD
    size_comm = comm.Get_size()
    rank = comm.Get_rank()
    
    nArguments = len(sys.argv)
    if nArguments != 4:
        text = """Three command line arguments are expected: \
                \n\t(1) ABQ inp file name including path \
                \n\t(2) Conectivity file name including path \
                \n\t(3) Output directory path"""
        raise RuntimeError(text)

    time_start = MPI.Wtime()    
    nameFileInpAbaqus = sys.argv[1]
    fileName_elementConnectivity = sys.argv[2]
    outputPath = sys.argv[3]
    
    #pathDir = os.path.splitdrive(nameFileInpAbaqus)[0] + os.path.split(os.path.splitdrive(nameFileInpAbaqus)[1])[0] + "/"
    pathDir = outputPath + "/"    
    # -------------------------------------------------------------------------------------------------------------------------------
    # -------------------------------------- Parse the ABQ inp file and obtain nodal coord information
    # -------------------------------------------------------------------------------------------------------------------------------
    listNodes, listCoordinatesX, listCoordinatesY, listCoordinatesZ = parseAbqInpFileForNodalCoords(nameFileInpAbaqus)
    listNodesRenumbered = [i+1 for i in range(0,len(listNodes))]
    # write the nodal information to a file
    out_path_nodesActual = pathDir + 'nodalCoords.dat'
    with open(out_path_nodesActual, 'w') as file_out_nodesActual:
        file_out_nodesActual.write(str(len(listNodes)) + '\n')
        for iNode in range(0,len(listNodes)):
            file_out_nodesActual.write("{0:8d}{1:25.10f}{2:25.10f}{3:25.10f}\n".format(listNodes[iNode], \
                                                                                       listCoordinatesX[iNode], \
                                                                                       listCoordinatesY[iNode], \
                                                                                       listCoordinatesZ[iNode]))
    file_out_nodesActual.close()
    
    out_path_nodesRenumbered = pathDir + 'nodalCoords_Renumbered.dat'
    with open(out_path_nodesRenumbered, 'w') as file_out_nodesRenumbered:
        file_out_nodesRenumbered.write(str(len(listNodesRenumbered)) + '\n')
        for iNode in range(0,len(listNodesRenumbered)):
            file_out_nodesRenumbered.write("{0:8d}{1:25.10f}{2:25.10f}{3:25.10f}\n".format(listNodesRenumbered[iNode],
                                                                                           listCoordinatesX[iNode],
                                                                                           listCoordinatesY[iNode],
                                                                                           listCoordinatesZ[iNode]
                                                                                           ))
    file_out_nodesRenumbered.close()
    # -------------------------------------------------------------------------------------------------------------------------------
    # -------------------------------------- Read connectivity information
    # -------------------------------------------------------------------------------------------------------------------------------
    listConnectivity = readConnectivity(fileName_elementConnectivity)
    
    time_start_ = MPI.Wtime()    
    startIndicesPartition, endIndicesPartition = performLinearPartitioning(len(listConnectivity),size_comm)
    listConnectivityPerProcessor = listConnectivity[startIndicesPartition[rank]:endIndicesPartition[rank]+1]      
    listConnectivityRenumberedPerProcessor = renumberConnectivity(listConnectivityPerProcessor,listNodes)
    comm.barrier()
    if rank == 0: print('Time for renumbering : ', MPI.Wtime()-time_start_,flush=True)
            
    if rank == 0:
        listConnectivityRenumbered = []
        for iData in range(0,len(listConnectivityRenumberedPerProcessor)):
            listConnectivityRenumbered.append(listConnectivityRenumberedPerProcessor[iData])
        if size_comm > 1:    
            for iRank in range(1,size_comm):
                data = comm.recv(source=iRank)
                for iData in range(0,len(data)):
                    listConnectivityRenumbered.append(data[iData])      
        # ------- write the renumbered data to file 
        out_path_connectivityRenumbered = pathDir + 'elementConnectivity_Renumbered.dat'
        with open(out_path_connectivityRenumbered, 'w') as file_out_connectivityRenumbered:
            for iElement in range(0,len(listConnectivityRenumbered)):
                file_out_connectivityRenumbered.write(str(listConnectivityRenumbered[iElement][0]) + '    ' + \
                                                      str(listConnectivityRenumbered[iElement][1]) + '    ' + \
                                                      str(listConnectivityRenumbered[iElement][2]) + '    ' + \
                                                      str(listConnectivityRenumbered[iElement][3]) + '    ' + \
                                                      str(listConnectivityRenumbered[iElement][4]) + '    ' + \
                                                      str(listConnectivityRenumbered[iElement][5]) + '    ' + \
                                                      str(listConnectivityRenumbered[iElement][6]) + '    ' + \
                                                      str(listConnectivityRenumbered[iElement][7]) + '\n')
        file_out_connectivityRenumbered.close()     
    else:    
        comm.send(listConnectivityRenumberedPerProcessor,dest=0)
        listConnectivityRenumbered = None

    comm.barrier()
    # send the renumbered connectivity to remaining processes        
    if size_comm > 1:
        listConnectivityRenumbered = comm.bcast(listConnectivityRenumbered, root = 0)
    
def renumberConnectivity(listConnectivity,listNodes):
    arrayConnectivityRenumbered = np.array(listConnectivity)
    for iNodeActual in range(0,len(listNodes)):
        arrayConnectivityRenumbered[arrayConnectivityRenumbered == listNodes[iNodeActual]] = iNodeActual + 1                        
    return arrayConnectivityRenumbered.tolist()    
    
if __name__ == "__main__":
    main()