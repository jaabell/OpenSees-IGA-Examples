
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
# nu12 = 0.055  # Poisson's ratio
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
               tagPlaneStress2, 1e1*E1, 1e1*E2, 0, 0, 0, rho)



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
        t = 2.43e-4
        thickness = [16*t]
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
# on bending strips

# NEED TO KNOW FIRST TWO BENDING STRIPS
for key in bendingStrips_data:
    connectedSurfacesTags=bendingStrips_data[key]['connectedSurfacesTags']
    if connectedSurfacesTags == [0, 1]:
        surfTag_01=bendingStrips_data[key]['surfTag']
    elif connectedSurfacesTags == [1, 0]:
        surfTag_10=bendingStrips_data[key]['surfTag']


ben0_tag = bendingStrips_data['bendingStrip_12']
ben1_tag = bendingStrips_data['bendingStrip_13']
surface_2_nodes=surfaces_data_loaded[surfTag_01]['nodeTags'][0,:] # Nodes on 00-01 (uv,uv)
surface_3_nodes=surfaces_data_loaded[surfTag_10]['nodeTags'][0,:] # Nodes on 00-01 (uv,uv)
nodesOnBase = np.concatenate([surface_0_nodes, surface_1_nodes,surface_2_nodes,surface_3_nodes])


for node in nodesOnBase:
    node=int(node)
    ops.fix(node,1,1,1)

# Should be equalDofDict[constrained]=retained
equalDofDict={}

for node in links.keys():
    equalDofDict[node] = reduceNodes(node,links)
    rN = equalDofDict[node]
    cN = int(node)

    # If it has a retianer node
    if rN != None:
        if ops.nodeCoord(rN) != ops.nodeCoord(cN):
                print('ERROR')
                exit()
        rN=int(rN)
        cN=int(cN)

        print('relation ', cN, rN)

        if cN not in nodesOnBase:
            ops.remove('sp', node, 1)
            ops.remove('sp', node, 2)
            ops.remove('sp', node, 3)
            ops.equalDOF(rN,cN,1,2,3)
    # else:
    #     if cN in nodesOnBase:
    #         print("Fixing node in base = ",cN)
    #         ops.fix(cN,1,1,1)




print('Finished equalDofing')



# Getting penalty factor
print("Getting penalty factor")

ops.wipeAnalysis()
# ops.constraints("Penalty", 1e13, 1e13)
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

# K_mean = np.max(K)
print("K_mean = ", K_mean)
penaltyFactor = 10**int(np.ceil(np.log10(K_mean))+8)
print("penaltyFactor = ", penaltyFactor)


# Starting gravity analysis
print("Starting gravity analysis")
ops.wipeAnalysis()  

# create TimeSeries
# ops.timeSeries("Constant", 1)
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
ops.numberer("RCM")
ops.system("UmfPack")
# ops.algorithm("NewtonLineSearch", 'type', 'Bisection')
# ops.algorithm("Broyden")
ops.algorithm("NewtonLineSearch")
# ops.algorithm("Newton")
# ops.algorithm("SecantNewton")
nSteps = 10
ops.integrator("LoadControl", 1.0/nSteps)
ops.analysis("Static")
ok = ops.analyze(nSteps)
if ok !=0:
    print('Analysis failed!')
    exit()
ops.loadConst('-time', 0.0) 
# ops.setTime(0.0) 

# exit()


print("Starting analysis")

# Clearing previously defined analysis parameters
ops.wipeAnalysis()  

# create TimeSeries for first row

# in X
# ops.timeSeries("Path", 1, '-time', *(t.tolist()), '-values', *(dispY.tolist()), '-prependZero')
ops.timeSeries("Path", 2, '-time', *(t.tolist()), '-values', *(dispY.tolist()))

# create a plain load pattern
ops.pattern("Plain", 2, 2)
# ops.pattern("UniformExcitation", 2, 2, '-disp', 2)





# Creating sp for nodes on base
for node in nodesOnBase:
    node=int(node)

    # ops.sp(nodeTag, dof, factor)  
    # ops.remove('sp', node, 2, 1)
    ops.remove('mp', node)
    ops.remove('sp', node, 2)
    ops.sp(node, 2, 1.0)   

# # Should be equalDofDict[constrained]=retained
# equalDofDict={}

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

#         ops.equalDOF(rN,cN,1,2,3)
# print('Finished equalDofing')


# Create test
# ops.test("EnergyIncr", 1.0e-7, 30, 1)
# ops.test("NormDispIncr", 1.0e-5, 50, 1)
ops.test("NormDispIncr", tol, 30, 1)
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
ops.constraints("Transformation") # Segmentation fault with this, maybe because of the equal dofs
# ops.constraints("Penalty", 1e14, 1e14)
# ops.constraints("Penalty", penaltyFactor, penaltyFactor)

# create integrator
# ops.integrator("Newmark", 0.5, 0.25)
# create integrator
rho_inf = 0.5       # Rho_inf ranges from 0 to 1
alpha_m = (2.0 - rho_inf) / (1.0 + rho_inf)
alpha_f = 1.0 / (1.0 + rho_inf)
gamma_t = 0.5 - alpha_f + alpha_m
beta_t = ((1.0 - alpha_f + alpha_m)**2) / 4.0
ops.integrator('GeneralizedAlpha', alpha_m, alpha_f)

# create SOE
# ops.system("FullGeneral")
ops.system("UmfPack")


# create analysis object
ops.analysis("Transient")



# print("\n\n\nPRINTING DOMAIN-----------------------")
# ops.printModel()
# print("\n\n\nDONE PRINTING DOMAIN-----------------------")


# Checking visualization to see if model is correct
fdef = 1e1
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
# container.render()

# Setting Rayleigh damping

print("Calculating Rayleigh damping")

Nmodes = 6 # or whatever
w2 = ops.eigen(Nmodes)

# Pick your modes and damping ratios
wi = w2[0]**0.5; zetai = 0.02 # 5% in mode 1
wj = w2[5]**0.5; zetaj = 0.018 # 2% in mode 6

A = np.array([[1/wi, wi],[1/wj, wj]])
b = np.array([zetai,zetaj])

x = np.linalg.solve(A,2*b)

#             M    KT  KI  Kn
ops.rayleigh(x[0],x[1],0.0,0.0)


# Create recorder
ops.recorder('Node', '-file', 'GNL8HZ.out', '-closeOnWrite', '-time', '-node', *[3100], '-dof', *[2], *['disp'])
# ops.recorder('Node', '-file', 'Node3100_Y_linear.out', '-closeOnWrite', '-time', '-node', *[3100], '-dof', *[2], *['disp'])

# algorithm = {1:'KrylovNewton', 2: 'SecantNewton' , 4: 'RaphsonNewton',5: 'PeriodicNewton', 6: 'BFGS', 7: 'Broyden', 8: 'NewtonLineSearch'}
# algorithm = {1:'SecantNewton', 2: 'NewtonLineSearch' , 3: 'BFGS',4: 'Broyden', 5: 'Newton'}
# algorithm = {1:'SecantNewton', 2: 'NewtonLineSearch', 3: "Newton"}
# algorithm = {1:'SecantNewton', 2: 'NewtonLineSearch', 3: "Newton", 4: 'KrylovNewton', 5: 'Broyden'}
# algorithm = {1:'SecantNewton', 2: 'NewtonLineSearch', 3: "Newton"}#, 4: 'KrylovNewton', 5: 'Broyden'}
algorithm = {1:'NewtonLineSearch', 2: 'SecantNewton', 3: "Newton"}#, 5: 'Broyden'}
# Performing the analysis

# nSteps = int(tMax/deltaT)
nSteps = len(t)
tFinal=t[-1]
current_time = 0
iCount = 0 
tCount=0
current_dt = float(deltaT_int(current_time))
surfName = 0
while current_time<tFinal:

    print("current_time = ", current_time )
    print("current_dt = ", current_dt)
    print("tFinal = ", tFinal)
    print(current_time/tFinal*100,"%")

    for i in algorithm:
        print("Using ", algorithm[i])

        if algorithm[i] == 'NewtonLineSearch':
            # ops.algorithm(algorithm[i], '-type', 'Bisection')
            ops.algorithm("NewtonLineSearch", '-type', 'Bisection')
        else:
            ops.algorithm(algorithm[i])

        ok = ops.analyze(1, current_dt)
        if (ok != 0):
            print("Analysis failed!!!")
            print("Changing algorithm")
            continue
        else:
            break

    if ok!=0 and iCount==0:
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
        print('ok = ', ok )
        iCount = 0
        tCount += 1
        # current_dt = deltaT
        # current_dt = deltaT_vector[tCount]
        current_time += current_dt
        current_dt = float(deltaT_int(current_time))


    print("Creating surface")
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

    fileName_stl = f'/home/rodrigo/FelipeElgueta/OpenSees-IGA-Examples/modeloAspaFondef/8HZ/figuresGNL8HZ/surf_{surfName}.stl'
    surfName+=1
    exchange.export_stl(container, fileName_stl,vertex_spacing=1)
