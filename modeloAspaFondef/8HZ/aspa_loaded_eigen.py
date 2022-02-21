

from scipy.sparse.linalg import eigsh
from geomdl.visualization import VisVTK, VisMPL
from geomdl import exchange_vtk
from geomdl import exchange
from edgeHandler import *
from scipy import interpolate
import numpy as np
import opensees as ops
from math import *

use_bendingStrip = True

deg2rad = pi / 180

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

# Prior to calibration data
GPa = 1e9
E1 = 45*GPa  # Young's modulus N/m^2 Pa
E2 = 10*GPa
nu12 = 0.3  # Poisson's ratio
nu21 = nu12*E2/E1  # Poisson's ratio
G12 = 5*GPa  # Shear modulus
rho = 2000.0  # *9.807 # kg/m^3
mm = 1.0 / 1000
t = 2.43e-4


# Real Data
# E1 = 56000e6  # Young's modulus N/m^2 Pa
# E2 = E1 * 0.500892857
# # nu12 = 0.055  # Poisson's ratio
# nu12 = 0.3  # Poisson's ratio
# nu21 = nu12*E2/E1  # Poisson's ratio
# G12 = 4118.3e6  # Shear modulus
# rho = 2900.0  # *9.807 # kg/m^3
# mm = 1.0 / 1000
# t = 2.43e-4


# Creating necessary plane stress
tagPlaneStress1 = 1
ops.nDMaterial("ElasticOrthotropicPlaneStress",
               tagPlaneStress1, E1, E2, nu12, nu21, G12, rho)

# Creating necessary plane stress
tagPlaneStress2 = 2
ops.nDMaterial("ElasticOrthotropicPlaneStress",
               tagPlaneStress2, 2e1*E1, 2e1*E2, nu12, nu21, 2e1*G12, rho)

# Material for "filling"
tagMat1 = 3
use_filling = False
ops.uniaxialMaterial('ENT', tagMat1, GPa*1e-2)

# [ 10.56004887  25.9524172   40.42553869  85.1434273  110.36874156
#  119.4921095  155.57179775 214.6428994  256.71883262 273.81615355] con filling

# w/2/pi:  [ 10.50330348  25.76180261  40.23545074  81.33397949 105.58803612
#  112.10607208 148.28041731 159.77466728 196.23314894 231.57888383] sin filling

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
        shellType = 'KLShell'

        if zone == 0 :
            θ = θ_0
            t = t_0
        elif zone == 1:
            θ = θ_1
            t = t_1
        elif zone == 2:
            θ = θ_2
            t = t_2
        elif zone == 3:
            θ = θ_2
            t = t_3
        elif zone == 4:
            θ = θ_2
            t = t_4
        elif zone == 5:
            θ = θ_2
            t = t_5

        matTags = [1]*len(θ)

        print(matTags)

        thickness = [t]*len(θ)

    else:
        θ = [0]
        matTags = [2] 
        thickness = [1*t]
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
    
    for dim in controlPts:  # Unweighting control pts
        for point in dim:
            point[0:3] /= point[3]

    

    # Creating a Patch in OpenSees
    ops.IGA("Patch", patchTag, nodeStartTag, p, q, nPtsU, nPtsV,
            "-type", shellType,
            "-nonLinearGeometry", 1,
            "-planeStressMatTags", *matTags,
            "-gFact", *gFact,
            "-theta", *θ,
            "-thickness", *thickness,
            "-uKnot", *uKnot, "-vKnot", *vKnot, "-controlPts", *controlPts.flatten())

    # Get the nodes on current patch
    nodesMap.append(
        np.arange(nodeStartTag, ops.getNodeTags()[-1] + 1).tolist())

    # Update patchTag, nodeStartTag and materialTag
    surface['ops_patchTag'] = patchTag
    lastElTag = ops.getEleTags()[-1]
    patchTag = lastElTag + 1 # next patch tag




# Info on moving base

# from matplotlib.pylab import *

# amplitude = 0.14 / 2 # 14 cm de carrito
# deltaT = 0.005
# Fs=1/deltaT

# ω_min = 0   # Hz 
# ω_max = 8   # Hz

# tMax = 2
# t=np.arange(0,tMax,deltaT)
# ω_in = np.linspace(ω_min,ω_max,len(t))
# phase_in = np.cumsum(ω_in/Fs)
# y=np.sin(2*np.pi*phase_in)


# t2=np.arange(tMax,2*tMax,deltaT)
# y2=np.sin(2*np.pi*ω_max*(t2-0*deltaT))
# t2-=deltaT


# t=np.concatenate([t,t2])
# dispY=np.concatenate([y,y2])*amplitude

# dispYpp=np.gradient(dispY)

# plot(t,dispY,'-b')
# plot(t,dispYpp, '-r')
# show()

from freq import *

amplitude = 0.14 / 2 # 14 cm de carrito
t_steady = 2
tMax = 4
ω_min = 0   # Hz 
ω_max = 8   # Hz
nPoints_accel=400
nPoints_steady=400
[deltaT_int,t,dispY] = generateOscillation(amplitude,t_steady,tMax,ω_min, ω_max,nPoints_accel,nPoints_steady)
dispYpp=np.gradient(dispY)
plot(t,dispYpp, '-r')
plot(t,dispY,'-b')
print(len(t))
plot(t,dispY,'og')
show()






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
nodesOnBase = np.concatenate([surface_0_nodes, surface_1_nodes,surface_2_nodes,surface_3_nodes])

# First and second row
# nodesOnBase = np.concatenate([surface_0_nodes, surface_1_nodes,surface_2_nodes,surface_3_nodes,surface_4_nodes,surface_5_nodes,surface_6_nodes,surface_7_nodes])



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
        if top_node_coords != bottom_node_coords and top_node_coords[2] >= 0.5:
            links_fill.append([top_nodes[i],bottom_nodes[i]])

eleTag = ops.getEleTags()[-1] + 1 # Get tag of last added material and start from there
i=0
for pair in links_fill:
    top_node = int(pair[0])
    bottom_node = int(pair[1])
    top_node_coords = ops.nodeCoord(top_node)
    bottom_node_coords = ops.nodeCoord(bottom_node)
    print("top_node Cords: ", top_node_coords)
    print("bottom_node Cords: ", bottom_node_coords,'\n')
    eleNodes = [top_node,bottom_node]

    A = 1e-2 # m2
    rho = 1e-3
    if use_filling:
        ops.element('Truss', eleTag, *eleNodes, A, tagMat1, '-rho', rho)
    eleTag += 1

    if top_node_coords[-1] == 1.5827112913356758:
        i+=1
        x_0 = top_node_coords[0]
        x_1 = bottom_node_coords[0]
        y_0 = top_node_coords[1]
        y_1 = bottom_node_coords[1]
        plot([x_0,x_1],[y_0,y_1],'-or')

print('i = ', i)

show()

print('LISTO')
# exit()

# create TimeSeries
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


print("Starting analysis")

# Create test
ops.test("NormDispIncr", 1.0e-8, 30, 1)

# create SOE
# ops.system("FullGeneral")
ops.system("UmfPack")

# create DOF number
# ops.numberer("Plain")
ops.numberer("RCM")

# create constraint handler
# ops.constraints("Plain")
# ops.constraints("Lagrange")
ops.constraints("Transformation")
# ops.constraints("Penalty",1e10,1e10)


ops.algorithm("Linear")
# ops.algorithm("Newton")
# ops.algorithm("SecantNewton")
# ops.algorithm("NewtonLineSearch", 'type', 'Bisection')
# ops.algorithm("NewtonLineSearch")
# ops.algorithm("ModifiedNewton")
# ops.algorithm("KrylovNewton")
# ops.algorithm("BFGS")
# ops.algorithm("Broyden")

# create integrator

nSteps = 1
ops.integrator("LoadControl", 1.0 / nSteps)

# create analysis object
ops.analysis("Static")

print("\n\n\nPRINTING DOMAIN-----------------------")
ops.printModel()
print("\n\n\nDONE PRINTING DOMAIN-----------------------")

# Checking visualization to see if model is correct
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

        weight = point[3]

        # Add deformation scaled by fdef
        point[0] = (cordX ) * weight
        point[1] = (cordY ) * weight
        point[2] = (cordZ ) * weight

    nPoints = surf.ctrlpts_size_u * surf.ctrlpts_size_v
    shape = np.array(compatibility.flip_ctrlpts2d(
        surf.ctrlpts2d[:])).shape
    controlPts = np.array(controlPts).reshape(shape)
    controlPts = np.array(compatibility.flip_ctrlpts2d(controlPts))

    surf.set_ctrlpts(controlPts.reshape(nPoints, 4).tolist(),
                     surf.ctrlpts_size_u, surf.ctrlpts_size_v)

# Visualizing deformed surface
container.sample_size = 40
for surf in container:
    surf.evaluate()

# Visualization configuration
container.vis = VisVTK.VisSurface(ctrlpts=False, legend=False)

# Render the aspa
container.render()


# Analyzing
for j in range(nSteps):
    print("=================================")
    print(f"Load step {j}")
    print("=================================")
    result = ops.analyze(1)
    if result != 0:
        break
        exit(-1)
    else:

        # Adding deformation to control points
        fdef = 2e2

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
        container.sample_size = 40
        for surf in container:
            surf.evaluate()

        # Visualization configuration
        container.vis = VisVTK.VisSurface(ctrlpts=False, legend=False)

        # Render the aspa
        container.render()




# for j in range(nSteps):
#     print("=================================")
#     print(f"Load step {j}")
#     print("=================================")
#     result = ops.analyze(1)
#     if result != 0:
#         break
#         exit(-1)
#     else:
#         # Adding deformation to control points
#         fdef = 1e2

#         for surface_name in surfaces_data_loaded:

#             surface=surfaces_data_loaded[surface_name]

#             nodeTags = surface["nodeTags"]
#             nPtsU = surface['nPtsU']
#             nPtsV = surface['nPtsV']

#             nodes=nodeTags.flatten().tolist()

#             controlPts = surface['ctrlpts']

#             controlPts = controlPts.reshape(
#                 [nPtsU * nPtsV, 4]).tolist()


#             for n in nodes:
#                 # Get node position in list
#                 indexN = nodes.index(n)
#                 point = controlPts[indexN]

#                 cordX=float(ops.nodeCoord(n, 1))
#                 cordY=float(ops.nodeCoord(n, 2))
#                 cordZ=float(ops.nodeCoord(n, 3))

#                 dispX=float(ops.nodeDisp(n, 1))
#                 dispY=float(ops.nodeDisp(n, 2))
#                 dispZ=float(ops.nodeDisp(n, 3))

#                 weight = point[3]

#                 # Add deformation scaled by fdef
#                 point[0] = (cordX + fdef * dispX) * weight
#                 point[1] = (cordY + fdef * dispY) * weight
#                 point[2] = (cordZ + fdef * dispZ) * weight

#             nPoints = nPtsU * nPtsV
#             shape = np.array(controlPts).shape
#             controlPts = np.array(controlPts).reshape(shape)
#             # controlPts = np.array(compatibility.flip_ctrlpts2d(controlPts))

#             surf.set_ctrlpts(controlPts.reshape(nPoints, 4).tolist(),
#                              nPtsU, nPtsV)

#         # Visualizing deformed surface
#         container.sample_size = 40
#         for surf in container:
#             surf.evaluate()

#         # Visualization configuration
#         container.vis = VisVTK.VisSurface(ctrlpts=False, legend=False)

#         # Render the aspa
#         container.render()

#         print("\nNext load step\n")


# Eigenvalues

print("Starting Eigen")

# ops.wipeAnalysis()

# # create SOE
# ops.system("FullGeneral")

# # create DOF number
# ops.numberer("Plain")

# # create constraint handler
# ops.constraints("Plain")

# ops.system("BandSPD")
# ops.integrator("Newmark", 0.5, 0.25)
# ops.system("FullGeneral")

nodes = ops.getNodeTags()

Nnodes = len(nodes)
Neigenvalues = 10  # arpack can only compute N-1 eigvals


w2s = ops.eigen(Neigenvalues)
# w2s = ops.eigen('-solver','-fullGenLapack',Neigenvalues)
# w2s = ops.eigen('-standard','-symmBandLapack',Neigenvalues)


order = np.argsort(w2s)
w2s = np.array(w2s, dtype=np.float64)[order]
w=np.sqrt(w2s)
print("w: ", w)



for i, w2 in enumerate(w2s):
    w = sqrt(abs(w2))
    f = w / 2 / pi
    T = 1 / f
    print(f"{i} {w2} {w} {f} {T} ")


phi = np.zeros((3 * Nnodes, Neigenvalues))

for i, n in enumerate(nodes):
    for j in range(Neigenvalues):
        phi[3 * i:3 * i + 3, j] = ops.nodeEigenvector(n, j + 1)

print(f"ϕ = {phi}")

w = np.sqrt(abs(w2s))
print("w/2/pi: ", w/2/pi)
for j in range(Neigenvalues):
    w = sqrt(abs(w2s[j]))
    f = w / 2 / pi
    print("=================================")
    print(f"Eigenvalue {j}")
    print(f"f = {f} Hz")
    print("=================================")

    # Adding deformation to control points
    fdef = 1e0
    for i in range(len(container)):
        surf = container[i]
        nodes = nodesMap[i]
        controlPts = np.array(compatibility.flip_ctrlpts2d(surf.ctrlpts2d[:])).reshape([surf.ctrlpts_size_v * surf.ctrlpts_size_u, 4]).tolist()
        for n in nodes:
            # Get node position in list
            indexN = nodes.index(n)
            point = controlPts[indexN]

            eigenVector=phi[:,order[j]]

            # Add deformation scaled by fdef
            weight = point[3]
            # point[0] = (ops.nodeCoord(n)[0] + fdef * eigenVector[3*(n-1)] )* weight
            # point[1] = (ops.nodeCoord(n)[1] + fdef * eigenVector[3*(n-1)+1] )* weight
            # point[2] = (ops.nodeCoord(n)[2] + fdef * eigenVector[3*(n-1)+2] )* weight
            # print("n: ", n)
            # print("j: ", j)
            # print("ops.nodeEigenvector(n,j): ", ops.nodeEigenvector(n,j+1))

            # From portwood digital
            point[0] = (ops.nodeCoord(n)[0] + fdef * ops.nodeEigenvector(n,j+1,1) )* weight
            point[1] = (ops.nodeCoord(n)[1] + fdef * ops.nodeEigenvector(n,j+1,2) )* weight
            point[2] = (ops.nodeCoord(n)[2] + fdef * ops.nodeEigenvector(n,j+1,3) )* weight

        nPoints = surf.ctrlpts_size_u * surf.ctrlpts_size_v
        shape = np.array(compatibility.flip_ctrlpts2d(surf.ctrlpts2d[:])).shape
        controlPts = np.array(controlPts).reshape(shape)
        controlPts = np.array(compatibility.flip_ctrlpts2d(controlPts))

        surf.set_ctrlpts(controlPts.reshape(nPoints, 4).tolist(), surf.ctrlpts_size_u, surf.ctrlpts_size_v)

    # Visualizing deformed surface
    container.sample_size = 40
    for surf in container:
        surf.evaluate()

    # Render the aspa
    container.render()

    for i in range(len(container)):
        surf = container[i]
        nodes = nodesMap[i]
        controlPts = np.array(compatibility.flip_ctrlpts2d(surf.ctrlpts2d[:])).reshape([surf.ctrlpts_size_v * surf.ctrlpts_size_u, 4]).tolist()
        for n in nodes:
            # Get node position in list
            indexN = nodes.index(n)
            point = controlPts[indexN]

            # Add deformation scaled by fdef
            weight = point[3]
            point[0] = (ops.nodeCoord(n)[0] ) * weight
            point[1] = (ops.nodeCoord(n)[1] ) * weight
            point[2] = (ops.nodeCoord(n)[2] ) * weight

        nPoints = surf.ctrlpts_size_u * surf.ctrlpts_size_v
        shape = np.array(compatibility.flip_ctrlpts2d(surf.ctrlpts2d[:])).shape
        controlPts = np.array(controlPts).reshape(shape)
        controlPts = np.array(compatibility.flip_ctrlpts2d(controlPts))

        surf.set_ctrlpts(controlPts.reshape(nPoints, 4).tolist(), surf.ctrlpts_size_u, surf.ctrlpts_size_v)

    print("\nNext Eigenvalue\n")