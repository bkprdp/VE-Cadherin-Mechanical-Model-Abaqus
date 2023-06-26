from odbAccess import *
from abaqusConstants import *
from collections import Counter
import math
from odbMaterial import *
from odbSection import *
import time as tt

from itertools import groupby

odb = openOdb(path='filename.odb')
myAssembly = odb.rootAssembly
mySteps = odb.steps

nodeDataDict = {}
i = 0
numNodes = 0

for instanceName in myAssembly.instances.keys():
    nodesData = myAssembly.instances[instanceName].nodes
    numNodes = numNodes + len(nodesData)
    elesData = myAssembly.instances[instanceName].elements
    for currNode in nodesData:
        coordX =  float(currNode.coordinates[0])
        coordY =  float(currNode.coordinates[1])
        coordZ =  float(currNode.coordinates[2])
        instName = currNode.instanceName
        nodeNum = currNode.label
        
        roundcoordX = round(coordX,5)
        roundcoordY = round(coordY,5)
        roundcoordZ = round(coordZ,5)

        strcoordX = str(roundcoordX)
        strcoordY = str(roundcoordY)
        strcoordZ = str(roundcoordZ)
        
        n_x = 0
        n_y = 0
        n_z = 0

        for i in range(len(strcoordX)):
            if strcoordX[i] == ".":
                n_x = i + 3
                break
        
        for i in range(len(strcoordY)):
            if strcoordY[i] == ".":
                n_y = i + 3
                break

        for i in range(len(strcoordZ)):
            if strcoordZ[i] == ".":
                n_z = i + 3
                break
        
        key = str([strcoordX[0:n_x],strcoordY[0:n_y],strcoordZ[0:n_z]])
        
        if key in nodeDataDict.keys():
            nodeDataDict[key].append([instName,nodeNum])
        else:
            nodeDataDict[key] = [[instName,nodeNum]]

biCellularNodesData = []
triCellularNodesData = []
triCellularNodesData1 = []

fileBiCellular = "bicellular_instance_nodes_3D_Mesh_1"
fileTriCellular = "tricellular_instance_nodes_3D_Mesh_1"

fbi = open(fileBiCellular,"w")
ftri = open(fileTriCellular,"w")

for nodeCoord in nodeDataDict.keys():
    if len(nodeDataDict[nodeCoord]) == 2 :
        biCellularNodesData.append(nodeDataDict[nodeCoord])
        fbi.write(str(nodeDataDict[nodeCoord]))
        fbi.write("\n")
    elif len(nodeDataDict[nodeCoord]) == 3 :
        triCellularNodesData.append(nodeDataDict[nodeCoord])
        ftri.write(str(nodeDataDict[nodeCoord]))
        ftri.write("\n")
        

print(len(biCellularNodesData))
print(len(triCellularNodesData))

fbi.close()
ftri.close()