from scipy.sparse.linalg import eigsh
from geomdl.visualization import VisVTK, VisMPL
from geomdl import exchange_vtk
from geomdl import exchange
from edgeHandler import *
from scipy import interpolate
import numpy as np
import opensees as ops
from math import *


# Model creation

use_bendingStrip = True

deg2rad = pi / 180

GPa = 1e9
MPa = 1e6

# materials

# Ply orientations according to zone
θ_0 = deg2rad*np.array([45.0, -45.0, 0.0, 90.0, 0.0, 90.0, 0.0, 90.0, 45.0, -45.0, 0.0, 90.0, 0.0, 90.0, 45.0, -45.0, 45.0, -45.0, 45.0, -45.0, 45.0, -45.0])
θ_1 = deg2rad*np.array([45.0, -45.0, 0.0, 90.0, 0.0, 90.0, 0.0, 90.0, 45.0, -45.0, 0.0, 90.0, 45.0, -45.0, 45.0, -45.0, 45.0, -45.0, 45.0, -45.0])
θ_2 = deg2rad*np.array([45.0, -45.0, 0.0, 90.0, 0.0, 90.0, 0.0, 90.0, 45.0, -45, -45.0, 45.0, -45.0, 45.0])

# Thickness according to zone
t_0 = 3.864e-4
t_1 = 2.575e-4
t_2 = 3.0714e-4
t_3 = 2.679e-4
t_4 = 2.536e-4
t_5 = 2.429e-4

# Elasticity according to zone, 61600 MPa inicial
E1_0 = 67975.15*MPa
E1_1 = 61600*MPa
E1_2 = 63224.28*MPa
E1_3 = 61600*MPa
E1_4 = 36242.57*MPa
E1_5 = 61600*MPa
# For loading dictionaries using pickle

from six.moves import cPickle as pickle #for performance

def load_dict(filename_):
    with open(filename_, 'rb') as f:
        ret_di = pickle.load(f)
    return ret_di


surfaces_data_loaded = load_dict('surfaces_data.pkl')

surfList=[]

for surface_name in surfaces_data_loaded:
    surface=surfaces_data_loaded[surface_name]
    # for key in surface:
    #     print (key, surface[key])
    # print('\n')



    nPtsU = surface['nPtsU']
    nPtsV = surface['nPtsV']
    uKnot = surface['uKnot']
    vKnot = surface['vKnot']
    p = surface["p"]
    q = surface["q"]


    nodeStartTag= surface['nodeStartTag']
    nodeEndTag= surface['nodeEndTag']

    nodeTags=np.arange(nodeStartTag,nodeEndTag+1)
    nodeTags=nodeTags.reshape([nPtsU,nPtsV],order = 'F')

    surface["nodeTags"] = nodeTags

    # print('nodeTags = \n', nodeTags, '\n')

    # controlPts = surface['ctrlpts'].tolist()
    controlPts = surface['ctrlpts'].reshape([nPtsU * nPtsV, 4]).tolist()

    # surf.set_ctrlpts(controlPts.reshape(nPoints, 4).tolist(), surf.ctrlpts_size_u, surf.ctrlpts_size_v)


    # Create a NURBS surface instance
    surf = NURBS.Surface()

    # Set degrees
    surf.degree_u = p
    surf.degree_v = q

    # Set control points
    surf.set_ctrlpts(controlPts, nPtsU, nPtsV)

    nPoints = surf.ctrlpts_size_u * surf.ctrlpts_size_v
    shape = np.array(compatibility.flip_ctrlpts2d(surf.ctrlpts2d[:])).shape
    controlPts = np.array(controlPts).reshape(shape)
    controlPts = np.array(compatibility.flip_ctrlpts2d(controlPts))
    surf.set_ctrlpts(controlPts.reshape(nPoints, 4).tolist(), surf.ctrlpts_size_u, surf.ctrlpts_size_v)

    # Set knot vectors
    # surf.knotvector_u = knotvector.generate(surf.degree_u, surf.ctrlpts_size_u)
    # surf.knotvector_v = knotvector.generate(surf.degree_v, surf.ctrlpts_size_v)

    surf.knotvector_u = uKnot
    surf.knotvector_v = vKnot

    surfType = surface['type']

    # if surfType == 'Patch':
    #     print(surface["zone"]) 

    if not use_bendingStrip:
        if surfType != 'Patch':
            continue

    surfList.append(surf)

container = multi.SurfaceContainer(surfList)
container.sample_size = 20

for surf in container:
    surf.evaluate()


# Visualization configuration
container.vis = VisVTK.VisSurface(ctrlpts=False)

# Render the aspa
container.render()





# CREATING OPENSEES MODEL

ops.wipe()
ops.model('basic', '-ndm', 3, '-ndf', 3)



Xt = 350.0*MPa
Xc = 350.0*MPa
Yt = 300.0*MPa
Yc = 300.0*MPa
S = 90.0*MPa
c1 = 4.0e-6
c2 = 30.0
c3 = 2.0e-6
c4 = 0.8
c5 = 80.0
c6 = 0.0
c7 = 0.0
c8 = 0.0
c9 = 0.0
b = 1.0

# Creating necessary plane stress for bending strips
tagPlaneStress2 = 1
ops.nDMaterial("ElasticOrthotropicPlaneStress",
               tagPlaneStress2, 1e1*E1_0, 1e1*E1_0, 0, 0, 0, 1)

# Material for "filling"
tagMat1 = 2
use_filling = False
ops.uniaxialMaterial('ENT', tagMat1, GPa*1e-2)

# Materials for patches
E_s = [E1_0,E1_1,E1_2,E1_3,E1_4,E1_5]
for (tag,E) in zip([3,4,5,6,7,8],E_s):
    E1 = E * 0.725
    E2 = E1 * 0.500892857
    nu12 = 0.055  # Poisson's ratio
    nu21 = nu12*E2/E1  # Poisson's ratio
    G12 = E1/(61600*MPa) * 4488*MPa   # Shear modulus
    rho = 3181.27*11.157/13.83742381290443 # kg/m^3, desired weight of 11.157 kg
    print('rho = ', rho)

    vonPaepeParams = [E1, E2, nu12, nu21, G12, rho, 
            Xt, Xc, Yt, Yc, S, c1, c2, c3, c4, c5, c6, c7, c8, c9, b]
    ops.nDMaterial("VonPapaDamage", tag, *vonPaepeParams)


gFact = [0, 0, 0]

patchTag = 1
nodesMap = []



# Creating patches
for surface_name in surfaces_data_loaded:
    surface=surfaces_data_loaded[surface_name]

    surfType = surface["type"]
    if surfType == 'Patch':

        # Checking zone
        zone = surface['zone']
        print('zone = ', zone)
        shellType = 'KLShell'

        if zone == 0 :
            θ = θ_0
            t = t_0
            tagPlaneStress = 3
        elif zone == 1:
            θ = θ_1
            t = t_1
            tagPlaneStress = 4
        elif zone == 2:
            θ = θ_2
            t = t_2
            tagPlaneStress = 5
        elif zone == 3:
            θ = θ_2
            t = t_3
            tagPlaneStress = 6
        elif zone == 4:
            θ = θ_2
            t = t_4
            tagPlaneStress = 7
        elif zone == 5:
            θ = θ_2
            t = t_5
            tagPlaneStress = 8

        # Creating necessary material tags
        matTags = [tagPlaneStress]*len(θ)

        print(matTags)

        thickness = [t]*len(θ)

        surface['nLayers'] = len(thickness)

    else:
        θ = [0]
        matTags = [tagPlaneStress2] 
        thickness = [10*t]
        shellType = "KLShell_BendingStrip"

    if not use_bendingStrip:
        if surfType != 'Patch':
            continue
    
    nPtsU = surface['nPtsU']
    nPtsV = surface['nPtsV']
    uKnot = surface['uKnot']
    vKnot = surface['vKnot']
    p = surface["p"]
    q = surface["q"]

    controlPts = surface['ctrlpts']

    nodeStartTag= surface['nodeStartTag']
    nodeEndTag = surface['nodeEndTag']
    
    for dim in controlPts:  # Unweighting control pts
        for point in dim:
            point[0:3] /= point[3]

    

    print("Start igaPatch")
    print("nodeStartTag = ", nodeStartTag)
    print("nodeEndTag = ", nodeEndTag)
    # Creating a Patch in OpenSees
    ops.IGA("Patch", patchTag, nodeStartTag, p, q, nPtsU, nPtsV,
            "-type", shellType,
            "-nonLinearGeometry", 1,
            "-planeStressMatTags", *matTags,
            "-gFact", *gFact,
            "-theta", *θ,
            "-thickness", *thickness,
            "-uKnot", *uKnot, "-vKnot", *vKnot, "-controlPts", *controlPts.flatten())
    print("Done igaPatch")


    # Get the nodes on current patch
    nodesMap.append(
        np.arange(nodeStartTag, ops.getNodeTags()[-1] + 1).tolist())

    # Update patchTag, nodeStartTag and materialTag
    surface['ops_patchTag'] = patchTag
    lastElTag = ops.getEleTags()[-1]
    surface['elementTags'] = np.arange(patchTag+1, lastElTag+1 ,1) # including lastElTag
    patchTag = lastElTag + 1 # next patch tag

    print('patchTag = ', patchTag)

    print("Finished surface[i]!\n")



print("Finished making model")

from matplotlib.pylab import *
from scipy.interpolate import interp1d

# 4_43 Hz case
# Getting data from input file
# input_path = '/home/rodrigo/FelipeElgueta/OpenSees-IGA-Examples/modeloAspaFondef/experimental/iga_inputs/case01_190715_13_20_05_y_base.txt'
input_path = '/home/rodrigo/FelipeElgueta/OpenSees-IGA-Examples/modeloAspaFondef/experimental/iga_inputs/case00_190609_11_21_47_ypp.txt'
Ypp_input = np.loadtxt(input_path,skiprows=1)
# nData_input = 3003
nData_input = len(Ypp_input)
dT_input = 0.00333
t_input = np.arange(0,nData_input*dT_input,dT_input)
Ypp_input_int = interp1d(t_input,Ypp_input)
Fs = 250
dT=1.0/Fs
start_time = 0
t = np.arange(start_time,t_input[-1],dT)
Ypp_4Hz = Ypp_input_int(t)

deltaT_vector = np.diff(t)
deltaT_vector = np.append(deltaT_vector,[deltaT_vector[-1]])
deltaT_int = interp1d(t,deltaT_vector)


# 7 Hz case
# Getting data from input file
# input_path = '/home/rodrigo/FelipeElgueta/OpenSees-IGA-Examples/modeloAspaFondef/experimental/iga_inputs/case01_190715_13_20_05_y_base.txt'
input_path = '/home/rodrigo/FelipeElgueta/OpenSees-IGA-Examples/modeloAspaFondef/experimental/iga_inputs/case01_190613_16_08_36_ypp.txt'
Ypp_input = np.loadtxt(input_path,skiprows=1)
# nData_input = 3003
nData_input = len(Ypp_input)
dT_input = 0.00333
t_input = np.arange(0,nData_input*dT_input,dT_input)
Ypp_input_int = interp1d(t_input,Ypp_input)
Fs = 250
dT=1.0/Fs
start_time = 0
t = np.arange(start_time,t_input[-1],dT)
Ypp_7Hz = Ypp_input_int(t)

deltaT_vector = np.diff(t)
deltaT_vector = np.append(deltaT_vector,[deltaT_vector[-1]])
deltaT_int = interp1d(t,deltaT_vector)


print("\n\n\nPRINTING DOMAIN-----------------------")
ops.printModel()
print("\n\n\nDONE PRINTING DOMAIN-----------------------")


# Looking for equalDofs

bendingStrips_data = load_dict('bendingStrip_dict.pkl')

links={}


for bendingStrip_name in bendingStrips_data:
    bendingStrip=bendingStrips_data[bendingStrip_name]

    print("bendingStrip_name = ", bendingStrip_name)

    # for key in bendingStrip:
    #     print (key, bendingStrip[key])

    retainedSurfaceTag = bendingStrip['connectedSurfacesTags'][0]
    constrainedSurfaceTag = bendingStrip['connectedSurfacesTags'][1]

  

    print('retainedSurfaceTag = ', retainedSurfaceTag)
    print("constrainedSurfaceTag = ", constrainedSurfaceTag)



    retainedEdge = bendingStrip['connectedEdges'][0]
    constrainedEdge = bendingStrip['connectedEdges'][1]

    # print('\n\n\n')
    if retainedSurfaceTag>constrainedSurfaceTag:
        retainedSurfaceTag, constrainedSurfaceTag = constrainedSurfaceTag, retainedSurfaceTag
        # retainedEdge, constrainedEdge = constrainedEdge, retainedEdge
        print("Swapping")
        swap = True
    else:
        swap = False

    retainedSurfaceName = f"surf_{retainedSurfaceTag}"
    constrainedSurfaceName = f"surf_{constrainedSurfaceTag}"

    retainedSurfaceNodes = surfaces_data_loaded[retainedSurfaceName]['nodeTags']
    constrainedSurfaceNodes = surfaces_data_loaded[constrainedSurfaceName]['nodeTags']


    print("retainedEdge = ", retainedEdge)
    print("constrainedEdge = ", constrainedEdge)


    # Getting nodes on retained edge:

    pointA = retainedEdge[0]
    pointB = retainedEdge[1]

    if pointA == "10" and pointB == "11":
        edgePoints_retained = retainedSurfaceNodes[-1,:]
        nextToEdgePoints_retained = retainedSurfaceNodes[-2,:]

    elif pointA == "00" and pointB == "10":
        edgePoints_retained = retainedSurfaceNodes[:,0]
        nextToEdgePoints_retained = retainedSurfaceNodes[:,1]

    elif pointA == "01" and pointB == "11":
        edgePoints_retained = retainedSurfaceNodes[:,-1]
        nextToEdgePoints_retained = retainedSurfaceNodes[:,-2]

    elif pointA == "00" and pointB == "01":
        edgePoints_retained = retainedSurfaceNodes[0,:]
        nextToEdgePoints_retained = retainedSurfaceNodes[1,:]


    # Getting nodes on constrained edge:

    pointA = constrainedEdge[0]
    pointB = constrainedEdge[1]

    if pointA == "10" and pointB == "11":
        edgePoints_constrained = constrainedSurfaceNodes[-1,:]
        nextToEdgePoints_constrained = constrainedSurfaceNodes[-2,:]

    elif pointA == "00" and pointB == "10":
        edgePoints_constrained = constrainedSurfaceNodes[:,0]
        nextToEdgePoints_constrained = constrainedSurfaceNodes[:,1]

    elif pointA == "01" and pointB == "11":
        edgePoints_constrained = constrainedSurfaceNodes[:,-1]
        nextToEdgePoints_constrained = constrainedSurfaceNodes[:,-2]

    elif pointA == "00" and pointB == "01":
        edgePoints_constrained = constrainedSurfaceNodes[0,:]
        nextToEdgePoints_constrained = constrainedSurfaceNodes[1,:]

    



    # Now creating bending strips equalDofs



    surface_name=f"surf_{bendingStrip_name.split('_')[1]}"

    bendingStripNodes = surfaces_data_loaded[surface_name]['nodeTags']

    if (pointA == "10" and pointB == "11") or (pointA == "00" and pointB == "01"):

        firstRow=bendingStripNodes[0,:]
        secondRow=bendingStripNodes[1,:]
        thirdRow=bendingStripNodes[2,:]

    else:

        firstRow=bendingStripNodes[:,0]
        secondRow=bendingStripNodes[:,1]
        thirdRow=bendingStripNodes[:,2]

    # secondRow va con retained Nodes
    # first row va con nextToEdgePoints_retained
    # third row va con nextToEdgePoints_constrained

    # print(len(firstRow),len(nextToEdgePoints_retained))
    # print(len(secondRow),len(edgePoints_retained))
    # print(len(thirdRow),len(nextToEdgePoints_constrained))

    firstRow_retained=nextToEdgePoints_retained
    secondRow_retained=edgePoints_retained
    thirdRow_retained=nextToEdgePoints_constrained

    if (retainedEdge == ['01', '11']) and (constrainedEdge == ['00','10']) and False: 

        print(" \n\n CASO RARO \n\n ")

        edgePoints_retained = retainedSurfaceNodes[:,0]
        nextToEdgePoints_retained = retainedSurfaceNodes[:,1]
        edgePoints_constrained = constrainedSurfaceNodes[:,-1]
        nextToEdgePoints_constrained = constrainedSurfaceNodes[:,-2]

        print("edgePoints_retained = ", edgePoints_retained, len(edgePoints_retained))
        print("edgePoints_constrained = ", edgePoints_constrained, len(edgePoints_constrained))

        firstRow_retained=nextToEdgePoints_retained
        secondRow_retained=edgePoints_retained
        thirdRow_retained=nextToEdgePoints_constrained


        # firstRow=bendingStripNodes[2,:]
        # secondRow=bendingStripNodes[1,:]
        # thirdRow=bendingStripNodes[0,:]

        firstRow=bendingStripNodes[:,2]
        secondRow=bendingStripNodes[:,1]
        thirdRow=bendingStripNodes[:,0]




    print("edgePoints_retained = ", edgePoints_retained, len(edgePoints_retained))
    print("edgePoints_constrained = ", edgePoints_constrained, len(edgePoints_constrained))

    for i in range(len(edgePoints_retained)):
        retainedNode=int(edgePoints_retained[i])
        constrainedNode=int(edgePoints_constrained[i])

        if ops.nodeCoord(retainedNode) != ops.nodeCoord(constrainedNode):
            print(ops.nodeCoord(retainedNode))
            print(ops.nodeCoord(constrainedNode))
            print('ERROR')
            print('\n')
            exit()

        # print(ops.nodeCoord(retainedNode))
        # print(ops.nodeCoord(constrainedNode))

        if constrainedNode in links:
            links[constrainedNode].append(retainedNode)
        else:
            links[constrainedNode] = [retainedNode] 




    if use_bendingStrip:

        print("\n\n\nComparing first row \n\n\n")
        for i in range(len(firstRow)):
            retainedNode=int(firstRow_retained[i])
            constrainedNode=int(firstRow[i])
            if swap == True:
                print("Swapping happened")
                # constrainedNode=int(thirdRow[i])

            if ops.nodeCoord(retainedNode) != ops.nodeCoord(constrainedNode):
                print(ops.nodeCoord(retainedNode))
                print(ops.nodeCoord(constrainedNode))
                print('ERROR')
                print('\n')
                exit()

            if constrainedNode in links:
                links[constrainedNode].append(retainedNode)
            else:
                links[constrainedNode] = [retainedNode]  

        print("\n\n\nComparing second row \n\n\n")
        for i in range(len(secondRow)):
            retainedNode=int(secondRow_retained[i])
            constrainedNode=int(secondRow[i])

            if ops.nodeCoord(retainedNode) != ops.nodeCoord(constrainedNode):
                print(ops.nodeCoord(retainedNode))
                print(ops.nodeCoord(constrainedNode))
                print('ERROR')
                print('\n')
                exit()

            if constrainedNode in links:
                links[constrainedNode].append(retainedNode)
            else:
                links[constrainedNode] = [retainedNode] 

        print("\n\n\nComparing third row \n\n\n")
        for i in range(len(thirdRow)):
            retainedNode=int(thirdRow_retained[i])
            constrainedNode=int(thirdRow[i])

            if ops.nodeCoord(retainedNode) != ops.nodeCoord(constrainedNode):
                print(ops.nodeCoord(retainedNode))
                print(ops.nodeCoord(constrainedNode))
                print('ERROR')
                print('\n')
                exit()

            if constrainedNode in links:
                links[constrainedNode].append(retainedNode)
            else:
                links[constrainedNode] = [retainedNode] 

        print("bendingStripNodes = \n",bendingStripNodes)


    # exit()
    print('\n')


# Checking if the paired nodes are correct
for key in links:
    # print(key,links[key])
    constrainedNode = key
    retainedNodes=links[key]
    for rN in retainedNodes:
        if ops.nodeCoord(rN) != ops.nodeCoord(constrainedNode):
            print(ops.nodeCoord(rN))
            print(ops.nodeCoord(constrainedNode))
            print('ERROR')
            exit()


# Adding null link in case the node has no links
for node in ops.getNodeTags():
    if node not in links:
        links[node]=[]




def reduceNodes(cN,links):
    # cN is the constrained node
    # links[cN] are the retainer nodes
    if links[cN]==[]:
        # equalDofDict[constrained] = None
        # print('end!',cN,links[cN])
        return None
    else:
        masterNode = None
        tempNodes = []
        for node in links[cN]:
            masterNode = reduceNodes(node,links)
            if masterNode != None:
                tempNodes.append(masterNode)
            else:
                masterNode = node
                tempNodes.append(masterNode)
        if len(tempNodes)>1:
            # print("constrainedNode = ",node, ops.nodeCoord(node))
            for node_tmp in tempNodes:
                if node_tmp != node:
                    masterNode = node_tmp
  
        return masterNode

# creating sp constraints on the nodes

# Nodes on the base
surface_0_nodes=surfaces_data_loaded['surf_0']['nodeTags'][0,:] # Nodes on 00-01 (uv,uv)
surface_1_nodes=surfaces_data_loaded['surf_1']['nodeTags'][0,:] # Nodes on 00-01 (uv,uv)


# Nodes on second row
surface_4_nodes=surfaces_data_loaded['surf_0']['nodeTags'][1,:] # Nodes on 00-01 (uv,uv)
surface_5_nodes=surfaces_data_loaded['surf_1']['nodeTags'][1,:] # Nodes on 00-01 (uv,uv)
# on bending strips

# NEED TO KNOW FIRST TWO BENDING STRIPS
for key in bendingStrips_data:
    connectedSurfacesTags=bendingStrips_data[key]['connectedSurfacesTags']
    if connectedSurfacesTags == [0, 1]:
        surfTag_01=bendingStrips_data[key]['surfTag']
    elif connectedSurfacesTags == [1, 0]:
        surfTag_10=bendingStrips_data[key]['surfTag']

surface_2_nodes=surfaces_data_loaded[surfTag_01]['nodeTags'][0,:] # Nodes on 00-01 (uv,uv)
surface_3_nodes=surfaces_data_loaded[surfTag_10]['nodeTags'][0,:] # Nodes on 00-01 (uv,uv)

# Second row
surface_6_nodes=surfaces_data_loaded[surfTag_01]['nodeTags'][1,:] # Nodes on 00-01 (uv,uv)
surface_7_nodes=surfaces_data_loaded[surfTag_10]['nodeTags'][1,:] # Nodes on 00-01 (uv,uv)


# Only first row
# nodesOnBase = np.concatenate([surface_0_nodes, surface_1_nodes,surface_2_nodes,surface_3_nodes])

# First and second row
nodesOnBase = np.concatenate([surface_0_nodes, surface_1_nodes,surface_2_nodes,surface_3_nodes,surface_4_nodes,surface_5_nodes,surface_6_nodes,surface_7_nodes])

def setEqualDofsAndBase():

    # Should be equalDofDict[constrained]=retained
    equalDofDict={}

    for node in links.keys():
        equalDofDict[node] = reduceNodes(node,links)
        rN = equalDofDict[node]
        cN = int(node)

        # If it has a retainer node
        if rN != None:
            if ops.nodeCoord(rN) != ops.nodeCoord(cN):
                    print('ERROR')
                    exit()
            rN=int(rN)
            cN=int(cN)

            # print('relation ', cN, rN)
            print(cN , ' retained to ', rN)

            if links[rN] != []:
                print('Node already retained!')
                exit()

            if cN not in nodesOnBase and rN not in nodesOnBase:
                # ops.remove('sp', cN, 1)
                # ops.remove('sp', cN, 2)
                # ops.remove('sp', cN, 3)
                # ops.remove('sp', rN, 1)
                # ops.remove('sp', rN, 2)
                # ops.remove('sp', rN, 3)
                ops.equalDOF(rN,cN,1,2,3)
        # else:
        #     if cN in nodesOnBase:
        #         print("Fixing node in base = ",cN)
        #         ops.fix(cN,1,1,1)

    for node in nodesOnBase:
        node=int(node)
        ops.fix(node,1,1,1)

    # exit()


    print('Finished equalDofing')

def createFilling():

    # Creating links for filling
    patchTags_list=[]
    i=0
    for surface_name in surfaces_data_loaded:
        surface=surfaces_data_loaded[surface_name]
        surfType = surface["type"]
        ops_patchTag = surface['ops_patchTag']
        if surfType == 'Patch':
            print('Patch ', ops_patchTag)
            patchTags_list.append(i)
            i+=1

    print(patchTags_list)
    patchTags_array=np.array(patchTags_list).reshape([int(len(patchTags_list)/2),2])
    print(patchTags_array)
    links_fill=[]
    for pair in patchTags_array:
        top_tag=pair[0]
        bottom_tag=pair[1]

        top_nodes = surfList[top_tag]['nodeTags']

        # top_nodes=surfaces_data_loaded[f'surf_{top_tag}']['nodeTags'][1:,:-1].flatten() # Nodes on 00-01 (uv,uv)
        # bottom_nodes=surfaces_data_loaded[f'surf_{bottom_tag}']['nodeTags'][1::,::-1].flatten() # Nodes on 00-01 (uv,uv)

        top_nodes=surfaces_data_loaded[f'surf_{top_tag}']['nodeTags'][1:,:].flatten() # Nodes on 00-01 (uv,uv)
        bottom_nodes=surfaces_data_loaded[f'surf_{bottom_tag}']['nodeTags'][1:,::-1].flatten() # Nodes on 00-01 (uv,uv)

        for i in range(len(top_nodes)):
            top_node = int(top_nodes[i])
            bottom_node = int(bottom_nodes[i])
            top_node_coords = ops.nodeCoord(top_node)
            bottom_node_coords = ops.nodeCoord(bottom_node)

            # 0.5 because thats when the filling starts, have to change that
            if top_node_coords != bottom_node_coords and top_node_coords[2] >= 0.0:
                links_fill.append([top_nodes[i],bottom_nodes[i]])

    eleTag = ops.getEleTags()[-1] + 1 # Get tag of last added element and start from there
    i=0
    for pair in links_fill:
        top_node = int(pair[0])
        bottom_node = int(pair[1])
        top_node_coords = ops.nodeCoord(top_node)
        bottom_node_coords = ops.nodeCoord(bottom_node)
        print("top_node Cords: ", top_node_coords)
        print("bottom_node Cords: ", bottom_node_coords,'\n')
        eleNodes = [top_node,bottom_node]

        A = 1e-3 # m2
        if use_filling:
            ops.element('Truss', eleTag, *eleNodes, A, tagMat1)

            # Just showing an example
            if top_node_coords[-1] == 1.5827112913356758:
                i+=1
                x_0 = top_node_coords[0]
                x_1 = bottom_node_coords[0]
                y_0 = top_node_coords[1]
                y_1 = bottom_node_coords[1]
                plot([x_0,x_1],[y_0,y_1],'-or')
        eleTag += 1


    print('i = ', i)

    show()



# Finished creating model


def getTol():
    # Getting penalty factor
    print("Getting penalty factor")

    ops.wipeAnalysis()
    ops.system('FullGeneral')
    ops.analysis('Transient')
     
    # Getting stiffness matrix
    ops.integrator('GimmeMCK',0.0,0.0,1.0)
    ops.analyze(1,0.0)
    # Number of equations in the model
    N = ops.systemSize() # Has to be done after analyze
    # Tolerance for norm disp (Taken from portwood digital)
    tol = 1e-5*np.sqrt(N)
    print("tol = ", tol)
    K = ops.printA('-ret')
    K = np.array(K)
    K.shape = (N,N)
    # K_mean = np.mean(K)
    K_mean = np.mean(np.diag(K))
    # K_mean = np.max(np.diag(K))

    # K_mean = np.max(K)
    print("K_mean = ", K_mean)
    penaltyFactor = 10**int(np.ceil(np.log10(K_mean))+8)
    print("penaltyFactor = ", penaltyFactor)

    return tol




def gravityAnalysis(nSteps):
    # Starting gravity analysis
    print("Starting gravity analysis")

    ops.wipeAnalysis()  

    ops.remove('timeSeries', 1)
    ops.remove('loadPattern', 1)

    ops.remove('timeSeries', 2)
    ops.remove('loadPattern', 2)

    for node in nodesOnBase:
        node=int(node)
        ops.remove('sp', node, 1)
        ops.remove('sp', node, 2)
        ops.remove('sp', node, 3)
        ops.fix(node,1,1,1)

    # Create timeSeries
    ops.timeSeries("Linear", 1)

    # create a plain load pattern
    ops.pattern("Plain", 1, 1)

    g = -9.807
    weight = [0, g, 0]

    # Loading patches
    for surface_name in surfaces_data_loaded:
        surface=surfaces_data_loaded[surface_name]
        
        patchType = surface['type']

        # Only load patches, not bendingStrips
        if patchType == 'Patch':
            patchTag = int(surface['ops_patchTag'])
            ops.eleLoad("-ele", patchTag, "-type", "-SelfWeight", *weight)

    # ops.test("NormDispIncr", 1.0e-5, 40, 1)
    ops.test("NormDispIncr", tol, 20, 1)
    # ops.constraints("Penalty", 1e13, 1e13)
    # ops.constraints("Plain")
    ops.constraints("Transformation")
    # ops.constraints("Penalty", penaltyFactor, penaltyFactor)
    # ops.constraints("Lagrange")
    ops.numberer("RCM")
    ops.system("UmfPack")
    # ops.system("SparseGEN")
    # ops.system("BandSPD")
    # ops.system("ProfileSPD")
    ops.algorithm("NewtonLineSearch", 'type', 'Bisection')

    # ops.algorithm("Linear")

    ops.integrator("LoadControl", 1.0/nSteps)

    ops.analysis("Static")

    ok = ops.analyze(nSteps)
    if ok !=0:
        print('Analysis failed!')
        exit()

    ops.loadConst('-time', 0.0)         


def rayleighDamping():
    # Setting Rayleigh damping

    print("Calculating Rayleigh damping")

    ops.wipeAnalysis()
    ops.system("UmfPack")
    # ops.system("ProfileSPD")
    ops.numberer("RCM")
    ops.constraints("Transformation")
    # ops.constraints("Penalty", penaltyFactor, penaltyFactor)
    ops.algorithm("Linear")


    Nmodes = 6 # or whatever
    ws = ops.eigen(Nmodes)
    w=np.sqrt(ws)
    f = w / 2 / np.pi
    print("f: ", f)


    # Pick your modes and damping ratios
    # wi = ws[0]**0.5; zetai = 0.02 # 5% in mode 1
    # wj = ws[5]**0.5; zetaj = 0.018 # 2% in mode 6

    wi = ws[0]**0.5; zetai = 0.8/100.0 # 5% in mode 1
    wj = ws[1]**0.5; zetaj = 1.0/100.0 # 2% in mode 6


    A = np.array([[1/wi, wi],[1/wj, wj]])
    b = np.array([zetai,zetaj])

    x = np.linalg.solve(A,2*b)

    # ops.modalDamping(1.005310)

    #             M    KT  KI  Kn
    ops.rayleigh(x[0],x[1],0.0,0.0)

    # ops.setTime(0.0) 

    # exit()

    return f

print("Starting analysis")

def resetAnalysisAndAddPathSp(inputCase):
    # Clearing previously defined analysis parameters
    ops.wipeAnalysis()  
    # ops.remove('timeSeries',2)

    # create TimeSeries 
    if inputCase==0:
        ops.timeSeries("Path", 2, '-time', *(t.tolist()), '-values', *(Ypp_4Hz.tolist()))
    
    else:
        ops.timeSeries("Path", 2, '-time', *(t.tolist()), '-values', *(Ypp_7Hz.tolist()))

    # create a UniformExcitation load pattern
    g = 9.807
    ops.pattern("UniformExcitation", 2, 2,'-accel',2,'-fact', g)


def analysisParameters():

    # wipe analysis
    ops.wipeAnalysis()  

    # ops.test("EnergyIncr", 1.0e-7, 30, 1)
    # ops.test("NormDispIncr", 1.0e-5, 50, 1)
    ops.test("NormDispIncr", tol, 12, 1)
    # ops.test("NormDispIncr", 1.0e-10, 2, 1)

    # ops.algorithm("Linear")
    # ops.algorithm("Newton")
    # ops.algorithm("SecantNewton")
    # ops.algorithm("NewtonLineSearch", 'type', 'Bisection')
    # ops.algorithm("NewtonLineSearch")
    # ops.algorithm("ModifiedNewton")
    # ops.algorithm("KrylovNewton")
    # ops.algorithm("BFGS")
    # ops.algorithm("Broyden")

    # create DOF number
    # ops.numberer("Plain")
    ops.numberer("RCM")

    # create constraint handler
    # ops.constraints("Plain") # Cannot use this one, does not allow for sp constrains other than fix
    # ops.constraints("Lagrange")
    # ops.constraints("Penalty", 1e14, 1e14)
    ops.constraints("Transformation") # Segmentation fault with this, maybe because of the equal dofs
    # ops.constraints("Penalty", penaltyFactor, penaltyFactor)

    # create integrator
    rho_inf = 0.5       # Rho_inf ranges from 0 to 1
    alpha_m = (2.0 - rho_inf) / (1.0 + rho_inf)
    alpha_f = 1.0 / (1.0 + rho_inf)
    gamma_t = 0.5 - alpha_f + alpha_m
    beta_t = ((1.0 - alpha_f + alpha_m)**2) / 4.0

    ops.integrator("Newmark", 0.5, 0.25)
    # ops.integrator('GeneralizedAlpha', alpha_m, alpha_f)

    # create SOE
    # ops.system("FullGeneral")
    ops.system("UmfPack")
    # ops.system("ProfileSPD")

    # create analysis object
    ops.analysis("Transient")

def resetMaxStress():
    # Loading patches
    for surface_name in surfaces_data_loaded:
        surface=surfaces_data_loaded[surface_name]
        
        patchType = surface['type']

        # Only load patches, not bendingStrips
        if patchType == 'Patch':
            patchTag = int(surface['ops_patchTag'])
            ops.setParameter('-val', 0, "-ele", patchTag, "resetMaxStress")
    pass



def fatigueAnalysis(i_cycle,t_Final,t_steadyState):

    ops.setTime(0.0)

    if t_Final > t[-1]:
        print("Error, tFinal too large")
        exit()

    # Folder to store surfaces in
    createFolder(i_cycle)

    # Create recorder
    # ops.remove('recorders')

    # lastNodeTag = ops.getNodeTags()[-1]
    lastNodeTag = 2320
    lastNodeTag = 3100
    print('lastNode coordinates: ', ops.nodeCoord(lastNodeTag) )
    # ops.recorder('Node', '-file', 'Linear_fatigue.out', '-closeOnWrite', '-time', '-node', *[4155], '-dof', *[2], *['disp'])

    #Useful recorders
    # ops.recorder('Node', '-file', f'Linear_fatigue_{i_cycle}.out', '-closeOnWrite', '-time', '-node', *[3033], '-dof', *[2], *['disp'])
    # ops.recorder('Node', '-file', f'Linear_fatigue_tip_{i_cycle}.out', '-closeOnWrite', '-time', '-node', *[2320], '-dof', *[2], *['disp'])

    ops.recorder('Node', '-file', f'GNL_fatigue_{i_cycle}.out', '-closeOnWrite', '-time', '-node', *[3033], '-dof', *[2], *['accel'])
    ops.recorder('Node', '-file', f'GNL_fatigue_tip_{i_cycle}.out', '-closeOnWrite', '-time', '-node', *[2097], '-dof', *[2], *['accel'])
    print('recorder 2320 coordinates: ', ops.nodeCoord(2320) )
    print('recorder 3033 coordinates: ', ops.nodeCoord(3033) )
    # exit()

    # ops.recorder('Node', '-file', 'Linear_fatigue_tip.out', '-closeOnWrite', '-time', '-node', *[2320], '-dof', *[2], *['disp'])
    # ops.recorder('Node', '-file', 'Node3100_Y_linear.out', '-closeOnWrite', '-time', '-node', *[3100], '-dof', *[2], *['disp'])

    # algorithm = {1:'KrylovNewton', 2: 'SecantNewton' , 4: 'RaphsonNewton',5: 'PeriodicNewton', 6: 'BFGS', 7: 'Broyden', 8: 'NewtonLineSearch'}
    # algorithm = {1:'SecantNewton', 2: 'NewtonLineSearch' , 3: 'BFGS',4: 'Broyden', 5: 'Newton'}
    # algorithm = {1:'SecantNewton', 2: 'NewtonLineSearch', 3: "Newton"}
    # algorithm = {1:'SecantNewton', 2: 'NewtonLineSearch', 3: "Newton", 4: 'KrylovNewton', 5: 'Broyden'}
    # algorithm = {1:'SecantNewton', 2: 'NewtonLineSearch', 3: "Newton"}#, 4: 'KrylovNewton', 5: 'Broyden'}
    algorithm = {1:'NewtonLineSearch', 2: 'SecantNewton', 3: "Newton", 4: 'KrylovNewton', 5:'Broyden'}
    # algorithm = {1:'Linear'}

    # ops.analysis("Transient")
    # Performing the analysis

    # nSteps = int(tMax/deltaT)
    nSteps = len(t)
    current_time = 0
    iCount = 0 
    tCount=0
    # current_dt = float(deltaT_int(current_time))
    current_dt = dT
    surfName = 0
    while current_time<t_Final:
        print("current_time = ", current_time )
        print("current_dt = ", current_dt)
        print("tFinal = ", t_Final)
        print(current_time/t_Final*100,"%")

        for i in algorithm:
            ops.algorithm(algorithm[i])
            print("Using ", algorithm[i])

            ok = ops.analyze(1, current_dt)
            if (ok != 0):
                print("Analysis failed!!!")
                print("Changing algorithm")
                continue
            else:
                break

        if ok!=0 and iCount<4:
            print("Analysis failed!!!")
            print("No more algorithms")
            print("Repeating with decreased DT")
            current_dt /= 2
            iCount += 1
            continue

        elif ok!=0 and iCount>0:
            print("Analysis failed!!!")
            exit()

        else:
            # print('ok = ', ok )
            iCount = 0
            tCount += 1
            # current_dt = deltaT
            # current_dt = deltaT_vector[tCount]
            current_time += current_dt
            # current_dt = float(deltaT_int(current_time))
            current_dt = dT

        # Reset max stress if t < t_steadyState
        if current_time <= t_steadyState:
            print("Resetting max stresses")
            resetMaxStress()
            pass
        else:
            print("Recording stresses!")

        # print("Creating surface")
        # current_time += current_dt
        # Adding deformation to control points
        fdef = 1e0

        for i in range(len(container)):
            surf = container[i]
            nodes = nodesMap[i]
            controlPts = np.array(compatibility.flip_ctrlpts2d(surf.ctrlpts2d[:])).reshape(
                [surf.ctrlpts_size_v * surf.ctrlpts_size_u, 4]).tolist()
            for n in nodes:
                # Get node position in list
                indexN = nodes.index(n)
                point = controlPts[indexN]

                cordX=float(ops.nodeCoord(n, 1))
                cordY=float(ops.nodeCoord(n, 2))
                cordZ=float(ops.nodeCoord(n, 3))

                dispX=float(ops.nodeDisp(n, 1))
                dispY=float(ops.nodeDisp(n, 2))
                dispZ=float(ops.nodeDisp(n, 3))


                weight = point[3]

                # Add deformation scaled by fdef
                point[0] = (cordX + fdef * dispX) * weight
                point[1] = (cordY + fdef * dispY) * weight
                point[2] = (cordZ + fdef * dispZ) * weight

            nPoints = surf.ctrlpts_size_u * surf.ctrlpts_size_v
            shape = np.array(compatibility.flip_ctrlpts2d(
                surf.ctrlpts2d[:])).shape
            controlPts = np.array(controlPts).reshape(shape)
            controlPts = np.array(compatibility.flip_ctrlpts2d(controlPts))

            surf.set_ctrlpts(controlPts.reshape(nPoints, 4).tolist(),
                             surf.ctrlpts_size_u, surf.ctrlpts_size_v)

        # Visualizing deformed surface
        container.sample_size = 10
        for surf in container:
            surf.evaluate()

        # Visualization configuration
        container.vis = VisVTK.VisSurface(ctrlpts=False, legend=False)

        # Render the aspa
        # container.render()

        # Export as vtk
        # fileName_vtk = f'figures/surf_{j}.vtk'
        # exchange_vtk.export_polydata(container, fileName_vtk, tessellate=True)

        fileName_stl = f'/home/rodrigo/FelipeElgueta/OpenSees-IGA-Examples/modeloAspaFondef/4_43Hz_accel/figuresLinear_fatigue_{i_cycle}/surf_{surfName}.stl'
        surfName+=1
        # exchange.export_stl(container, fileName_stl,vertex_spacing=1)

    pass


def eigen(nModes):
    print("Calculating Eigen")

    ops.wipeAnalysis()
    ops.system("UmfPack")
    # ops.system("ProfileSPD")
    ops.numberer("RCM")
    ops.constraints("Transformation")
    # ops.constraints("Penalty", penaltyFactor, penaltyFactor)
    ops.algorithm("Linear")

    ws = ops.eigen(nModes)
    w=np.sqrt(ws)
    f = w / 2 / np.pi
    print("f: ", f)

    return f


def createFolder(i_cycle):

    import os

    parent_dir = '/home/rodrigo/FelipeElgueta/OpenSees-IGA-Examples/modeloAspaFondef/4_43Hz_accel/'

    # directory = f'figuresLinear_fatigue_{i_cycle}'
    directory = f'figuresLinear_fatigue_{i_cycle}'

    path = os.path.join(parent_dir,directory)

    try: 
        os.mkdir(path)
    except OSError as error:
        print(error)  


def getNjump(nJump_perc,max_incr_dam):

    # Just in case
    nJump = 1

    nJump_list = [] #don't know the final size

    # Loading patches
    for surface_name in surfaces_data_loaded:
        surface=surfaces_data_loaded[surface_name]
        
        patchType = surface['type']

        # Only load patches, not bendingStrips
        if patchType == 'Patch':
            # patchTag = int(surface['ops_patchTag'])
            elementTags = [int(tag) for tag in surface['elementTags']] 
            nGauss_patch = int((surface['p']+1)*(surface['q']+1))
            nLayers_patch = surface['nLayers']

            # looping through elements with damage model
            for (eleTag,eleNumber) in zip(elementTags,range(len(elementTags))):
                for gp in range(nGauss_patch): 
                    for layer in range(nLayers_patch):
                        nJump_ele = ops.eleResponse(eleTag, "material_and_layer_number", str(gp), str(layer), 'NJUMP', str(max_incr_dam))
                        # nJump_ele = ops.eleResponse(eleTag, "material_and_layer_number", str(gp), str(layer), 'calculateNJUMP', str(nJump_perc))

                        for nJump_type in nJump_ele:
                            nJump_list.append(nJump_type)

    # Creating histogram
    hist, bins = np.histogram(nJump_list, bins=int(1e6), density=True)

    cumulativeFreq = np.cumsum(np.diff(bins)*hist) # This is the integral of the distribution
    cumulativeFreq = np.insert(cumulativeFreq,0,[0.0])

    # index_njump_perc = (np.abs(cumulativeFreq-nJump_perc)).argmin()
    # print('index_njump_perc = ', index_njump_perc)
    # print("int(bins[index_njump_perc]) = ", int(bins[index_njump_perc]))

    for j in range(len(cumulativeFreq)):
        # print("cumulativeFreq[j] = ", cumulativeFreq[j])
        if cumulativeFreq[j] >= nJump_perc:
            index_njump_perc = j-1
            print("int(bins[index_njump_perc]) = ", int(bins[index_njump_perc]))
            break

    frequency_inter=interp1d(cumulativeFreq,bins)
    nJump = int(frequency_inter(nJump_perc)) # number of cycles to advance, interpolating

    if nJump == 0 :
        print("nCycles = 0 ")
        nJump = 1

    # plot(bins,cumulativeFreq)
    # plot(nJump,nJump_perc,'ob')
    # show()

    return nJump



def advanceDamageState(nJump):
    print("Advancing damage state in ", nJump, "cycles")
    # Loading patches
    for surface_name in surfaces_data_loaded:
        surface=surfaces_data_loaded[surface_name]
        
        patchType = surface['type']

        # Only load patches, not bendingStrips
        if patchType == 'Patch':
            patchTag = int(surface['ops_patchTag'])
            ops.setParameter('-val', int(nJump), "-ele", patchTag, "advanceDamageState") # Only calling to each patch
    pass

def writeOnFile(string):
    # Writing to file
    with open("modalResults_GNL_test.txt","a") as file:
        # Writing data to a file
        file.write(string_result)


# Setting model
setEqualDofsAndBase()
createFilling()
tol = getTol()

# Creating MPCO recorder
ops.recorder("mpco","aspa_GNL",
    "-N","displacement",
    "-E","section.fiber.stress",
    "-E","section.fiber.strain",
    "-E","section.fiber.damagestate",
    )

nCycles = 0
i_cycle = 0
nModes = 6

file_modal_Results = open("modalResults_GNL_test.txt",'a')

# First step 
nSteps = 10
f = eigen(nModes)

# Writing frequencies and cycles on file
string_result = f'{nCycles} '
for mode in f:
    string_result += f' {mode} '
string_result += '\n'

# file_modal_Results.write(string_result)
writeOnFile(string_result)


# First gravity analysis
gravityAnalysis(nSteps)

g = 9.807
bladeWeight = 0
ops.reactions()
for node in nodesOnBase:
    node=int(node)
    bladeWeight += ops.nodeReaction(node, 2)

bladeWeight /= g
print('bladeWeight [kg] = ', bladeWeight)


# Get cycles data from txt
cyclesData = np.loadtxt('cycles_data.txt', skiprows = 1)
nCycles_max = int(cyclesData[-1,-1])

nCases = len(cyclesData[:,0])
inputCase = 0 # 0 if 4.4 Hz, 1 if 7.5 Hz

i=0
nCycles_this = int(cyclesData[i,5])
while nCycles < nCycles_max:
    nCycles_this = int(cyclesData[i,5]) # number of cycles for this test
    print(f"Case {i}: {nCycles_this} cycles")
    while nCycles < nCycles_this:

        print(f"{nCycles} cycles")
        print(f"{nCycles/nCycles_this*100} %")

        # Do analysis with corresponding input_case
        resetAnalysisAndAddPathSp(inputCase) # 0 if 4.4 Hz, 1 if 7.5 Hz

        t_Final = 3.3
        t_steadyState = 3

        f = rayleighDamping()
        analysisParameters()
        fatigueAnalysis(i_cycle, t_Final, t_steadyState)

        # Once the cycles finish, Calculate njump vector, get nJump and advance damage state
        nJump_perc = 0.02
        max_incr_dam = 1e-2
        nJump = getNjump(nJump_perc, max_incr_dam)

        # Check if calculated nJump is bound to this test, reduce if that's the case
        if (nCycles+nJump)>=nCycles_this:
            nJump = nCycles_this - nCycles

        advanceDamageState(nJump) 

        nCycles += nJump # Update this with the number of calculated cycles
        i_cycle += 1 # Update this with 1 (number of current "jump")

        # Redo eigen analysis
        f = eigen(nModes)

        # Writing frequencies and cycles on file
        string_result = f'{nCycles} '
        for mode in f:
            string_result += f' {mode} '
        string_result += '\n'

        # file_modal_Results.write(string_result)
        writeOnFile(string_result)

        # Resetting model
        ops.reset()

        # Setting time 0
        ops.setTime(0.0)

        # Setting up next step
        gravityAnalysis(10)

    # Next test case
    i+=1
    if inputCase == 0:
        inputCase = 1
    else:
        inputCase = 0 




    

exit()











resetAnalysisAndAddPathSp()

while nCycles < nCycles_max:

    t_Final = 4.3
    t_steadyState = 4

    f = rayleighDamping()
    analysisParameters()
    fatigueAnalysis(i_cycle, t_Final, t_steadyState)

    # Once the cycles finish, Calculate njump vector, get nJump and advance damage state
    nJump_perc = 0.02
    max_incr_dam = 1e-2
    nJump = getNjump(nJump_perc, max_incr_dam)
    advanceDamageState(nJump) 

    nCycles += nJump # Update this with the number of calculated cycles
    i_cycle += 1 # Update this with 1 (number of current "jump")

    # Redo eigen analysis
    f = eigen(nModes)

    # Writing frequencies and cycles on file
    string_result = f'{nCycles} '
    for mode in f:
        string_result += f' {mode} '
    string_result += '\n'

    # file_modal_Results.write(string_result)
    writeOnFile(string_result)

    # Resetting model
    ops.reset()

    # Setting time 0
    ops.setTime(0.0)

    # Setting up next step
    gravityAnalysis(10)
    resetAnalysisAndAddPathSp()