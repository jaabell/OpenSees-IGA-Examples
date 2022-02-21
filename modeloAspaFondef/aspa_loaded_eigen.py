
from scipy.sparse.linalg import eigsh
from geomdl.visualization import VisVTK, VisMPL
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

# Real Data
# E1 = 56000e6  # Young's modulus N/m^2 Pa
E1 = 45000.0e6  # Young's modulus N/m^2 Pa
# E1 = 58300e6  # Young's modulus N/m^2 Pa
# E2 = E1 * 0.500892857
E2 = 10000.0e6
nu12 = 0.3  # Poisson's ratio
# nu12 = 0.055  # Poisson's ratio
nu21 = nu12*E2/E1  # Poisson's ratio
# G12 = 4488.0e6  # Shear modulus
G12 = 5000.0e6  # Shear modulus
# G12 = 4118.3e6  # Shear modulus
# rho = 3181.27  # *9.807 # kg/m^3
rho = 2000.0  # *9.807 # kg/m^3
# rho = 2900.0  # *9.807 # kg/m^3
mm = 1.0 / 1000



# Creating necessary plane stress
tagPlaneStress1 = 1
ops.nDMaterial("ElasticOrthotropicPlaneStress",
               tagPlaneStress1, E1, E2, nu12, nu21, G12, rho)

# Creating necessary plane stress
tagPlaneStress1 = 2
ops.nDMaterial("ElasticOrthotropicPlaneStress",
               tagPlaneStress1, 1e1*E1, 1e1*E2, 0, 0, 0, rho)



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
        matTags = [2] * len(θ)
        thickness = [t_0]*len(θ)
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
            # "-nonLinearGeometry", 0,
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

# Fixing base nodes
ops.fixZ(-0.229, 1, 1, 1)  # First row
ops.fixZ(-0.119, 1, 1, 1)  # Second row



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

    retainedEdge = bendingStrip['connectedEdges'][0]
    constrainedEdge = bendingStrip['connectedEdges'][1]


    retainedSurfaceName = f"surf_{retainedSurfaceTag}"
    constrainedSurfaceName = f"surf_{constrainedSurfaceTag}"

    retainedSurfaceNodes = surfaces_data_loaded[retainedSurfaceName]['nodeTags']
    constrainedSurfaceNodes = surfaces_data_loaded[constrainedSurfaceName]['nodeTags']



    print('retainedSurfaceTag = ', retainedSurfaceTag)
    print("constrainedSurfaceTag = ", constrainedSurfaceTag)

    print("retainedEdge = ", retainedEdge)
    print("constrainedEdge = ", constrainedEdge)

    # print('\n\n\n')


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

    if (retainedEdge == ['01', '11']) and (constrainedEdge == ['00','10']): 

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
        print(ops.nodeCoord(retainedNode))
        print(ops.nodeCoord(constrainedNode))
        print('\n')

        if constrainedNode in links:
            links[constrainedNode].append(retainedNode)
        else:
            links[constrainedNode] = [retainedNode] 




    if use_bendingStrip:

        print("\n\n\nComparing first row \n\n\n")
        for i in range(len(firstRow)):
            retainedNode=int(firstRow_retained[i])
            constrainedNode=int(firstRow[i])

            print(ops.nodeCoord(retainedNode))
            print(ops.nodeCoord(constrainedNode))
            print('\n')

            if constrainedNode in links:
                links[constrainedNode].append(retainedNode)
            else:
                links[constrainedNode] = [retainedNode]  

        print("\n\n\nComparing second row \n\n\n")
        for i in range(len(secondRow)):
            retainedNode=int(secondRow_retained[i])
            constrainedNode=int(secondRow[i])

            print(ops.nodeCoord(retainedNode))
            print(ops.nodeCoord(constrainedNode))
            print('\n')

            if constrainedNode in links:
                links[constrainedNode].append(retainedNode)
            else:
                links[constrainedNode] = [retainedNode] 

        print("\n\n\nComparing third row \n\n\n")
        for i in range(len(thirdRow)):
            retainedNode=int(thirdRow_retained[i])
            constrainedNode=int(thirdRow[i])

            print(ops.nodeCoord(retainedNode))
            print(ops.nodeCoord(constrainedNode))
            print('\n')

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
        return None
    else:
        masterNode = None
        for node in links[cN]:
            if reduceNodes(node,links) != []:
                masterNode = node
        return masterNode
            


# Should be equalDofDict[constrained]=retained
equalDofDict={}

for node in links.keys():
    equalDofDict[node] = reduceNodes(node,links)
    rN = equalDofDict[node]
    cN = node

    if rN != None:
        if ops.nodeCoord(rN) != ops.nodeCoord(cN):
                print('ERROR')
                exit()
        rN=int(rN)
        cN=int(cN)

        ops.equalDOF(rN,cN,1,2,3)

print('Finished equalDofing')


# from mayavi import mlab

# fig =mlab.figure()

# x = []
# y=[]
# z=[]

# for node in links.keys():
#     equalDofDict[node] = reduceNodes(node,links)
#     rN = equalDofDict[node]
#     cN = node

#     if rN != None:
#         if ops.nodeCoord(rN) != ops.nodeCoord(cN):
#                 print('ERROR')
#                 exit()
#         rN=int(rN)
#         cN=int(cN)

#         xCoord_rN = ops.nodeCoord(rN,1)
#         yCoord_rN = ops.nodeCoord(rN,2)
#         zCoord_rN = ops.nodeCoord(rN,3)

#         xCoord_cN = ops.nodeCoord(cN,1) + 1e0
#         yCoord_cN = ops.nodeCoord(cN,2) + 1e0
#         zCoord_cN = ops.nodeCoord(cN,3) + 1e0

#         x.append(xCoord_rN)
#         y.append(yCoord_rN)
#         z.append(zCoord_rN)

#         x.append(xCoord_cN)
#         y.append(yCoord_cN)
#         z.append(zCoord_cN)

#         mlab.points3d(xCoord_rN, yCoord_rN, zCoord_rN, scale_factor=0.002, figure = fig)
#         mlab.points3d(xCoord_cN, yCoord_cN, zCoord_cN, scale_factor=0.002, figure = fig)

#         mlab.plot3d([xCoord_rN, xCoord_cN], [yCoord_rN, yCoord_cN], [zCoord_rN, zCoord_cN], figure = fig, line_width = 0.002)




        # mlab.plot3d(xCoord_rN, yCoord_rN, zCoord_rN)
        # mlab.plot3d(xCoord_cN, yCoord_cN, zCoord_cN)


        # ops.equalDOF(rN,cN,1,2,3)

# mlab.points3d(x, y, z, s, colormap="RdYlBu", scale_factor=0.02,
#               scale_mode='none', mode='2dcross')


# n = 5000
# x = np.random.rand(n)
# y = np.random.rand(n)
# z = np.random.rand(n)
# s = np.sin(x)**2 + np.cos(y)

# mlab.points3d(x, y, z, scale_factor=0.002)

# input()

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
ops.constraints("Plain")
# ops.constraints("Lagrange")
# ops.constraints("Transformation")
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