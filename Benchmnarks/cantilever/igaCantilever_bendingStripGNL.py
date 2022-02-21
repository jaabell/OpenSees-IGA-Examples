#  IGA CANTILEVER PLATE with bending strip

# For timing
import datetime

import numpy as np
import opensees as ops
from math import *



# Geomgl utilities for visualization and surface manipulation
from geomdl import NURBS, compatibility, operations, multi, knotvector
from surfVisualize import *

La = 5.0  	#
Lb = 1.0  	#
mm = 1.0 / 1000.  # m

ops.wipe()
ops.model('basic', '-ndm', 3, '-ndf', 3)

# These are given in v,u
controlPts = np.array([
  [0.00000  , 0.0 , 0.0 , 1.0] ,
  [0.00000  , 0.5 , 0.0 , 1.0] ,
  [0.00000  , 1.0 , 0.0 , 1.0] ,
  [0.31250  , 0.0 , 0.0 , 1.0] ,
  [0.31250  , 0.5 , 0.0 , 1.0] ,
  [0.31250  , 1.0 , 0.0 , 1.0] ,
  [0.93750  , 0.0 , 0.0 , 1.0] ,
  [0.93750  , 0.5 , 0.0 , 1.0] ,
  [0.93750  , 1.0 , 0.0 , 1.0] ,
  [1.87500  , 0.0 , 0.0 , 1.0] ,
  [1.87500  , 0.5 , 0.0 , 1.0] ,
  [1.87500  , 1.0 , 0.0 , 1.0] ,
  [3.12500  , 0.0 , 0.0 , 1.0] ,
  [3.12500  , 0.5 , 0.0 , 1.0] ,
  [3.12500  , 1.0 , 0.0 , 1.0] ,
  [4.37500  , 0.0 , 0.0 , 1.0] ,
  [4.37500  , 0.5 , 0.0 , 1.0] ,
  [4.37500  , 1.0 , 0.0 , 1.0] ,
  [5.62500  , 0.0 , 0.0 , 1.0] ,
  [5.62500  , 0.5 , 0.0 , 1.0] ,
  [5.62500  , 1.0 , 0.0 , 1.0] ,
  [6.87500  , 0.0 , 0.0 , 1.0] ,
  [6.87500  , 0.5 , 0.0 , 1.0] ,
  [6.87500  , 1.0 , 0.0 , 1.0] ,
  [8.12500  , 0.0 , 0.0 , 1.0] ,
  [8.12500  , 0.5 , 0.0 , 1.0] ,
  [8.12500  , 1.0 , 0.0 , 1.0] ,
  [9.06250  , 0.0 , 0.0 , 1.0] ,
  [9.06250  , 0.5 , 0.0 , 1.0] ,
  [9.06250  , 1.0 , 0.0 , 1.0] ,
  [9.68750  , 0.0 , 0.0 , 1.0] ,
  [9.68750  , 0.5 , 0.0 , 1.0] ,
  [9.68750  , 1.0 , 0.0 , 1.0] ,
  [10.00000 , 0.0 , 0.0 , 1.0] ,
  [10.00000 , 0.5 , 0.0 , 1.0] ,
  [10.00000 , 1.0 , 0.0 , 1.0]
])

for point in controlPts:
    point[0]/=2



# Create a NURBS surface instance
surfL = NURBS.Surface()

# Set surface degrees
surfL.degree_u = 4
surfL.degree_v = 1

# Setting control points for surface
surfL.set_ctrlpts(controlPts.tolist(), 12, 3)

# Set knot vectors
surfL.knotvector_u = knotvector.generate(surfL.degree_u, surfL.ctrlpts_size_u)
surfL.knotvector_v = knotvector.generate(surfL.degree_v, surfL.ctrlpts_size_v)


surfR = operations.translate(surfL, [La, 0, 0])


# Creating bending strip

from edgeHandler import *

interfacePoints, nodesOnLeftPatch=edgeGetter(surfL,"10","11")

nodesOnRightPatch = edgeGetter(surfR,"00","01")[1]

bendingStrip=makeBendingStrip(nodesOnLeftPatch,interfacePoints,nodesOnRightPatch,direction='u')



# Creating container for multipatches

surfList = [surfL, surfR, bendingStrip]

container = multi.SurfaceContainer(surfList)

# Visualize surface

container.sample_size = 5
for surf in container:
    surf.evaluate()
# Visualization configuration
container.vis = VisMPL.VisSurface(ctrlpts=True, legend=False)
# container.vis.ctrlpts_offset=0.1

# Render the surface
evalcolor = ["green", "green", 'blue']
container.render(evalcolor=evalcolor)

first_time=datetime.datetime.now()


# nDMaterial ElasticIsotropic $nDtag_elastic $elasticidad_probeta
# $poisson_probeta
E1 = 1.2e6  # Young's modulus N/m^2
E2 = E1
nu12 = 0.0  # Poisson's ratio
nu21 = 0.0  # Poisson's ratio
G12 = E1/(2*(1+nu21)) # Shear modulus
rho = 8.0e3  # *9.807 # kg/m^3
t = 100.0 * mm

MPa = 1e6
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

vonPaepeParams = [E1, E2, nu12, nu21, G12, rho, 
            Xt, Xc, Yt, Yc, S, c1, c2, c3, c4, c5, c6, c7, c8, c9, b]


tagPlaneStress1 = 1
# ops.nDMaterial("ElasticOrthotropicPlaneStress", tagPlaneStress1, E1, E2, nu12, nu21, G12, rho)
ops.nDMaterial("VonPapaDamage", tagPlaneStress1, *vonPaepeParams)

tagPlaneStress2 = 2
# ops.nDMaterial("ElasticOrthotropicPlaneStress", tagPlaneStress2, E1, E2, nu12, nu21, G12, rho)
ops.nDMaterial("VonPapaDamage", tagPlaneStress2, *vonPaepeParams)

tagPlaneStress3 = 3
bStripFactor = 1e2
# ops.nDMaterial("ElasticOrthotropicPlaneStress", tagPlaneStress3, bStripFactor*E1, bStripFactor*E1 , 0, 0, 0, rho)
vonPaepeParams = [bStripFactor*E1, bStripFactor*E2, 0, 0, 0, 0.01, 
            Xt, Xc, Yt, Yc, S, c1, c2, c3, c4, c5, c6, c7, c8, c9, b]
ops.nDMaterial("VonPapaDamage", tagPlaneStress3, *vonPaepeParams)


deg2rad = pi / 180
Nlayers = 3

matTags = [1]*Nlayers
thickness = [t/Nlayers]*Nlayers
θ = [0 * deg2rad]*Nlayers

gFact = [0.0, 0.0, 0.0 * 9.807]


Nlayers = len(θ)

names = ["leftPatch", "rightPatch", "bendingStrip"]
shellType = 'KLShell'
patchTag = 1
nodeStartTag = 1

nodesMap=[]


nLayers_list=[2,2,2]

for i in range(len(container)):
    surf = container[i]
    name = names[i]

    Nlayers = nLayers_list[i]

    matTags = [1]*Nlayers
    thickness = [t/Nlayers]*Nlayers
    θ = [0 * deg2rad]*Nlayers

    # Identifying bending strip
    if name == "bendingStrip":
        shellType = "KLShell_BendingStrip"

    # Flipping control point to u,v format
    controlPts = surf.ctrlpts2d[:]
    controlPts = np.array(compatibility.flip_ctrlpts2d(controlPts))

    for dim in controlPts: # Unweighting control pts
        for point in dim:
            point[0:3]/=point[3]

    # Creating a Patch in OpenSees
    ops.IGA("Patch", patchTag, nodeStartTag, surf.degree_u, surf.degree_v, surf.ctrlpts_size_u, surf.ctrlpts_size_v,
            "-type", shellType,
            # "-nonLinearGeometry", 0,
            "-planeStressMatTags", *matTags,
            "-gFact", *gFact,
            "-theta", *θ,
            "-thickness", *thickness,
            "-uKnot", *surf.knotvector_u, "-vKnot", *surf.knotvector_v, "-controlPts", *controlPts.flatten())


    # Get the nodes on current patch 
    nodesMap.append(np.arange(nodeStartTag,ops.getNodeTags()[-1]+1).tolist())

    # Update patchTag, nodeStartTag and materialTag
    lastElTag = ops.getEleTags()[-1]
    lastNodeTag = ops.getNodeTags()[-1]
    matTags[0] += 1
    patchTag = lastElTag + 1
    nodeStartTag = lastNodeTag + 1



# Creating constraints
for n in ops.getNodeTags():
    n = int(n)
    if n in [1, 2, 13, 14, 25, 26]:
        ops.fix(n, 1, 1, 1)


# Creating equal dof's in patch interface
interface1=[12, 24, 36]
interface2=[37, 49, 61]

for i in range(len(interface1)):
    ret = interface1[i]
    const = interface2[i]
    print('ret = ', ops.nodeCoord(ret))
    print('const = ', ops.nodeCoord(const))
    ops.equalDOF(ret, const, 1, 2, 3)

# nodesOnPatches = [11, 12, 25, 23, 24, 38]
nodesOnPatches = [11, 12, 38, 23, 24, 50, 35, 36, 62]
nodesOnBendingStrip = [73, 74, 75, 76, 77, 78, 79, 80 , 81]

for i in range(len(nodesOnPatches)):
    const = int(nodesOnPatches[i])
    ret = int(nodesOnBendingStrip[i])
    print('ret = ', ops.nodeCoord(ret))
    print('const = ', ops.nodeCoord(const))
    ops.equalDOF(const, ret, 1, 2, 3)


print("\n\n\nPRINTING DOMAIN-----------------------")
ops.printModel()
print("\n\n\nDONE PRINTING DOMAIN-----------------------")




# ------------------------------
# Start of analysis generation
# ------------------------------

# create TimeSeries
ops.timeSeries("Linear", 1)

# create a plain load pattern
ops.pattern("Plain", 1, 1)


print("Loading nodes")

# Pz = 4.0*1e-12
Pz = 4.0
for n in [48, 60, 72]:
    # # Vertical load
    if n == 60:
      ops.load(n, 0, 0, -Pz / 2.0)
    else:
      ops.load(n, 0, 0, -Pz / 4.0)

    # # Axial load
    # ops.load(n, -Pz/2.0, 0, 0)

I = Lb * (sum(thickness)**3) / 12.0

print("Finished loading nodes")


print("Starting analysis")

# Create test
ops.test("NormDispIncr", 1.0e-8, 500, 1)


# create SOE
ops.system("UmfPack")

# create DOF number
# ops.numberer("Plain")
ops.numberer("RCM")

# create constraint handler
# ops.constraints("Plain")
ops.constraints("Transformation")
# ops.constraints("Penalty",1,1)


# ops.algorithm("Linear")
# ops.algorithm("Newton")
# ops.algorithm("SecantNewton")
ops.algorithm("NewtonLineSearch")
# ops.algorithm("ModifiedNewton")
# ops.algorithm("KrylovNewton")
# ops.algorithm("BFGS")
# ops.algorithm("Broyden")

# create integrator

# nSteps = 50
# ops.integrator("LoadControl", 1.0 / nSteps)

delta = -0.4
defMax = 6.7
nSteps = abs(int(defMax / delta))
# ops.integrator("LoadControl", 1.0 / nSteps)
ops.integrator("DisplacementControl", 48, 3, delta)


# create analysis object
ops.analysis("Static")

# Create MPCO recorder
ops.recorder("mpco","iga_cantilever",
    "-N","displacement",
    "-E","section.fiber.stress",
    "-E","section.fiber.strain",
    "-E","section.fiber.damagestate",
    )

print(f"DONE! ")

import matplotlib.pyplot as plt
data = np.zeros((nSteps + 1, 3))

# perform the analysis
for j in range(nSteps):
    print("=================================")
    print(f"Load step {j}")
    print("=================================")
    result = ops.analyze(1)
    if result != 0:
        break
        exit(-1)
    data[j + 1, 0] = abs(ops.nodeDisp(60, 3))
    data[j + 1, 1] = abs(ops.nodeDisp(60, 1))
    data[j + 1, 2] = abs(ops.getLoadFactor(1) * Pz)


later_time=datetime.datetime.now()

diference=later_time - first_time
deltaT = divmod(diference.total_seconds(), 60)
print('deltaT = ',deltaT)


# Adding deformation to control points
fdef = 1e0
for i in range(len(container)):
    surf=container[i]
    nodes=nodesMap[i]
    print("nodes: ", nodes)
    controlPts = np.array(compatibility.flip_ctrlpts2d(surf.ctrlpts2d[:])).reshape([surf.ctrlpts_size_v * surf.ctrlpts_size_u, 4]).tolist()
    for n in nodes:
        # Get node position in list
        indexN = nodes.index(n)
        point = controlPts[indexN]

        # Add deformation scaled by fdef
        weight=point[3]
        point[0] = (ops.nodeCoord(n)[0] + fdef * ops.nodeDisp(n)[0])*weight
        point[1] = (ops.nodeCoord(n)[1] + fdef * ops.nodeDisp(n)[1])*weight
        point[2] = (ops.nodeCoord(n)[2] + fdef * ops.nodeDisp(n)[2])*weight


    nPoints=surf.ctrlpts_size_u*surf.ctrlpts_size_v
    shape=np.array(compatibility.flip_ctrlpts2d(surf.ctrlpts2d[:])).shape
    controlPts=np.array(controlPts).reshape(shape)
    controlPts=np.array(compatibility.flip_ctrlpts2d(controlPts))

    surf.set_ctrlpts(controlPts.reshape(nPoints,4).tolist(), surf.ctrlpts_size_u, surf.ctrlpts_size_v)



container.sample_size = 10
for surf in container:
    surf.evaluate()
    
# Visualization configuration
container.vis = VisMPL.VisSurface(ctrlpts=True, legend=False, animate=True)
# container.vis = VisMPL.VisSurfWireframe(ctrlpts=True, legend=False)
# container.vis.ctrlpts_offset=0.1

# Render the surface
evalcolor = ["green", "green", 'black']
cpcolor=["green","green","black"]
container.render(evalcolor=evalcolor, cpcolor=cpcolor)

result = np.array(ops.nodeDisp(60, 3))
print("result: ", result)

F=Pz*np.arange(0,1.05,0.05)
U=[0,0.026,0.103,0.224,0.381,0.563,0.763,0.971,1.184,1.396,1.604,1.807,2.002,2.19,2.37,2.541,2.705,2.861,3.01,3.151,3.286]
W=[0,0.663,1.309,1.922,2.493,3.015,3.488,3.912,4.292,4.631,4.933,5.202,5.444,5.66,5.855,6.031,6.19,6.335,6.467,6.588,6.698]
plt.plot(data[:, 0], data[:, 2], 'xr',label="W disp (IGA)") #w
plt.plot(W, F, '-r',label="W disp (reference)") #w

plt.plot(data[:, 1], data[:, 2], 'xb',label="U disp (GNL)") #u
plt.plot(U, F, '-b',label="U disp (reference)") #u

plt.xlabel('Horizontal Displacement')
plt.ylabel('Horizontal Load')
plt.legend()
plt.show()

print("Done")


print("Done!")

ops.remove("recorders")




# PostProcessing

import scipy.io as sio
import h5py

# file = h5py.File('aspa.mpco', 'r')
file = h5py.File('/home/rodrigo/FelipeElgueta/OpenSees-IGA-Examples/Benchmnarks/cantilever/iga_cantilever.mpco', 'r')

use_damage = True


surfList = [surfL, surfR]



stages = [1]
maxTimes = 16
startTime = 0

for stage in stages:

    modelStage = f'MODEL_STAGE[{stage}]'

    patches = file[modelStage]['MODEL']['ELEMENTS']['1-UnknownMovableObject[0:0]'][:, 0] # These are the subdomains -> Iga Patches
    elements = file[modelStage]['MODEL']['ELEMENTS']['207-UnknownMovableObject[0:0]'][:, 0] # These are the IGA_KLShell elements

    nodes = file[modelStage]['MODEL']['ELEMENTS']['207-UnknownMovableObject[0:0]'][:, 1:]

    print('nodes = \n', nodes)

    numPatches = 1

    nodesOnPatches = []
    elementsOnPatches = []

    continueSearch = True
    while continueSearch:

      for i in range(1,len(nodes)):
        currentElement = elements[i-1]
        nextElement = elements[i]

        print("currentElement = ", currentElement)
        print("nextElement = ", nextElement)

        foundPatch = True

        if abs(nextElement - currentElement) > 1:
          continueSearch = True
          foundPatch = True
          numPatches += 1
          print('breaking in i = ', i)
          break_i = i
          nodesOnPatches.append(nodes[:i,:])
          elementsOnPatches.append(elements[:i])
          nodes = nodes[i:,:]
          elements = elements[i:]
          print('nodes = ', nodes)
          print('elements = ', elements)
          break
        else:
          foundPatch = False

      if not foundPatch:
        continueSearch = False
        nodesOnPatches.append(nodes)
        elementsOnPatches.append(elements)

      else:
        pass


    print('numPatches = ', numPatches)

    # Getting layup info
    layupData = file[modelStage]['MODEL']['SECTION_ASSIGNMENTS']
    layupData_keys = layupData.keys()
    numLayups = len(layupData)
    print('numLayups = ', numLayups)
    print('layupData_keys = ', layupData_keys)

    zLocs_list = []

    for key in layupData_keys:
        print('key = ', key)

        layupData_here = layupData[key]

        zLocs_layups = layupData_here['FIBER_DATA'][:,1] # Shape (numMaterial x 3 (algo, zLoc, thickness))

        print('zLocs_layups = ', zLocs_layups)

        zLocs_list.append(zLocs_layups)




    # Getting nodeLimits in order to know which nodes are from which patch
    nodeLimits = []
    for nodes in nodesOnPatches:
      nodeLimits.append([nodes[0,0],nodes[-1,-1]])




    nTimes = len(file[modelStage]['RESULTS']['ON_NODES']['DISPLACEMENT']['DATA'])

    nTimes = min(maxTimes,nTimes) # Just so the data is not too much

    for time in range(nTimes):

      step = f'STEP_{time + startTime}' 

      controlPts = file[modelStage]['MODEL']['NODES']['COORDINATES'][:, :] # All the control points / nodes in the model
      U = file[modelStage]['RESULTS']['ON_NODES']['DISPLACEMENT']['DATA'][step][:,:] # All the displacements in the model

      nResults = len(file[modelStage]['RESULTS']['ON_ELEMENTS']['section.fiber.stress']) # Doesn't have to be stress, just need the number of files

      # The files are different because they are not the same size (nLayers)
     
      dictResults = {}

      for nFile in range(nResults):
        # Getting current file tag

        actualFile = f'207-UnknownMovableObject[0:0:{nFile}]'

        multiplicity = file[modelStage]['RESULTS']['ON_ELEMENTS']['section.fiber.stress'][actualFile]['META']['MULTIPLICITY'][0,0]


        dictResults[multiplicity] = actualFile
    

      controlPts_patches = []
      sigma_patches = []
      strain_patches =[]
      damage_patches = []
      U_patches=[]
      Nlayers_patches=[]


      # Looping through the 'patches'
      prevNlayers = len(zLocs_list[0])
      for i in range(len(nodeLimits)):

        zLocs = zLocs_list[i]

        Nlayers = len(zLocs)
        Nlayers_patches.append(Nlayers)

        # This is [nlayers1, nlayers2, nlayers3,.... etc]
        sizesFile = dictResults.keys()

        # for Nlayers_file in sizesFile:
        #     print('Nlayers_file = ', Nlayers_file)  
        #     if Nlayers == Nlayers_file

   

        # Need the specific file on strains and stresses depending on the size
        actualFile = dictResults[Nlayers]
        print('actualFile = ', actualFile)


        sigma_global = file[modelStage]['RESULTS']['ON_ELEMENTS']['section.fiber.stress'][actualFile]['DATA'][step][:,:] # Shape (numElements x Nlayers*noGps*3)
        strain_global = file[modelStage]['RESULTS']['ON_ELEMENTS']['section.fiber.strain'][actualFile]['DATA'][step][:,:] # Shape (numElements x Nlayers*noGps*3)
        
        if use_damage:
            damage_global = file[modelStage]['RESULTS']['ON_ELEMENTS']['section.fiber.damagestate'][actualFile]['DATA'][step][:,:] # Shape (numElements x Nlayers*noGps*9)

        nodeStart, nodeEnd = nodeLimits[i]
        rangeNodes = np.arange(nodeStart,nodeEnd+1,1) - 1
        grid = np.ix_(rangeNodes,range(3))

        controlPts_here = controlPts[grid]
        U_here = U[grid]
        shapeU_here = np.shape(U_here)
        U_here = np.reshape(U_here,shapeU_here[0]*shapeU_here[1])


        elementsThisPatch = elementsOnPatches[i]

        controlPts_patches.append(controlPts_here)
        U_patches.append(U_here)

        # print('elementsThisPatch = ', elementsThisPatch)

        if Nlayers != prevNlayers or i==0:
            # I'm on a different file
            newFile = True
        else:
            # Same file
            newFile = False



        # Have to do a new grid
        if newFile:
          # I'm in the first patch
          numEleThisPatch = len(elementsThisPatch)
          id_elements = np.arange(0,numEleThisPatch,1)

        else:
          elementsPrevPatch = elementsOnPatches[i-1]
          numEleThisPatch = len(elementsThisPatch)
          numElePrevPatch = len(elementsPrevPatch)
          id_elements = np.arange(numElePrevPatch, numElePrevPatch+numEleThisPatch,1)

        # print('id_elements = ', id_elements)
        print('len(id_elements) = ', len(id_elements))


        nData = len(sigma_global[0])
        grid = np.ix_(id_elements,range(nData))

        if use_damage:
            nData_damage = len(damage_global[0])
            grid_damage = np.ix_(id_elements,range(nData_damage))

        # continue





        # Getting stresses and strains from elements in each layer
        sigma = sigma_global[grid]
        strain = strain_global[grid]
        if use_damage:
            damage = damage_global[grid_damage]

        # print('sigma = ', sigma)


        shapeSigma=np.shape(sigma) # numEle x 3*nGauss
        shapeStrain=np.shape(strain)
        if use_damage:
            shapeDamage=np.shape(damage)

        numEle = shapeSigma[0] 

        # print('numEle = ', numEle)

        # Ranges for stress and strain
        range0 = np.arange(0,shapeSigma[1],3)
        range1 = range0+1
        range2 = range0+2

        # Ranges for damage
        if use_damage:
            range0_damage = np.arange(0,shapeDamage[1],9)
            range1_damage = range0_damage+1
            range2_damage = range0_damage+2
            range3_damage = range0_damage+3
            range4_damage = range0_damage+4
            range5_damage = range0_damage+5
            range6_damage = range0_damage+6
            range7_damage = range0_damage+7
            range8_damage = range0_damage+8

        # print('range0 = ', range0)

        # Getting quantities in each direction for all layers
        sigmaXX = sigma[:,np.ix_(range0)]
        sigmaYY = sigma[:,np.ix_(range1)]
        sigmaXY = sigma[:,np.ix_(range2)]
        strainXX = strain[:,np.ix_(range0)]
        strainYY = strain[:,np.ix_(range1)]
        strainXY = strain[:,np.ix_(range2)]

        # Getting damage variables (9 in total)
        if use_damage:
            damage_0 = damage[:,np.ix_(range0_damage)]
            damage_1 = damage[:,np.ix_(range1_damage)]
            damage_2 = damage[:,np.ix_(range2_damage)]
            damage_3 = damage[:,np.ix_(range3_damage)]
            damage_4 = damage[:,np.ix_(range4_damage)]
            damage_5 = damage[:,np.ix_(range5_damage)]
            damage_6 = damage[:,np.ix_(range6_damage)]
            damage_7 = damage[:,np.ix_(range7_damage)]
            damage_8 = damage[:,np.ix_(range8_damage)]


        # print('sigmaXX = ', sigmaXX)

        nGauss_3 = int(shapeSigma[1]/Nlayers)
        if use_damage:
            nGauss_9 = int(shapeDamage[1]/Nlayers)

        # print('nGauss_3 = ', nGauss_3)

        sigma = np.zeros([Nlayers,numEle,nGauss_3])
        strain = np.zeros([Nlayers,numEle,nGauss_3])
        if use_damage:
            damage = np.zeros([Nlayers,numEle,nGauss_9])



        for iLayer in range(Nlayers):

          rangeLayer = np.arange(iLayer,int(nGauss_3/3)*Nlayers,Nlayers)
          if use_damage:
            rangeLayer_damage = np.arange(iLayer,int(nGauss_9/9)*Nlayers,Nlayers)
          # if Nlayers == 1:
          #   rangeLayer = np.arange(iLayer,int(nGauss_3/3),Nlayers)
          # else:
          #   rangeLayer = np.arange(iLayer,nGauss_3,Nlayers)

          # print('rangeLayer = ', rangeLayer)

          # Getting stresses and strains corresponding to actual layer 
          sigmaXX_layer = sigmaXX[:,:,np.ix_(rangeLayer)]
          sigmaYY_layer = sigmaYY[:,:,np.ix_(rangeLayer)]
          sigmaXY_layer = sigmaXY[:,:,np.ix_(rangeLayer)]


          # print("sigmaXX_layer = ", sigmaXX_layer)
          # print("sigmaYY_layer = ", sigmaYY_layer)
          # print("sigmaXY_layer = ", sigmaXY_layer)
          # exit()


          strainXX_layer = strainXX[:,:,np.ix_(rangeLayer)]
          strainYY_layer = strainYY[:,:,np.ix_(rangeLayer)]
          strainXY_layer = strainXY[:,:,np.ix_(rangeLayer)]

          # Getting damage variables corresponding to actual layer
          if use_damage:
            damage0_layer = damage_0[:,:,np.ix_(rangeLayer_damage)]
            damage1_layer = damage_1[:,:,np.ix_(rangeLayer_damage)]
            damage2_layer = damage_2[:,:,np.ix_(rangeLayer_damage)]
            damage3_layer = damage_3[:,:,np.ix_(rangeLayer_damage)]
            damage4_layer = damage_4[:,:,np.ix_(rangeLayer_damage)]
            damage5_layer = damage_5[:,:,np.ix_(rangeLayer_damage)]
            damage6_layer = damage_6[:,:,np.ix_(rangeLayer_damage)]
            damage7_layer = damage_7[:,:,np.ix_(rangeLayer_damage)]
            damage8_layer = damage_8[:,:,np.ix_(rangeLayer_damage)]




          shapeSigma = np.shape(sigmaXX_layer)
          shapeStrain = np.shape(strainXX_layer)
          if use_damage:
            shapeDamage = np.shape(damage0_layer)

          # print('shapeSigma = ', shapeSigma)

          sigmaXX_layer = np.reshape(sigmaXX_layer,[1,shapeSigma[0],shapeSigma[-1]])
          sigmaYY_layer = np.reshape(sigmaYY_layer,[1,shapeSigma[0],shapeSigma[-1]])
          sigmaXY_layer = np.reshape(sigmaXY_layer,[1,shapeSigma[0],shapeSigma[-1]])

          strainXX_layer = np.reshape(strainXX_layer,[1,shapeStrain[0],shapeStrain[-1]])
          strainYY_layer = np.reshape(strainYY_layer,[1,shapeStrain[0],shapeStrain[-1]])
          strainXY_layer = np.reshape(strainXY_layer,[1,shapeStrain[0],shapeStrain[-1]])

          if use_damage:
            damage0_layer = np.reshape(damage0_layer,[1,shapeDamage[0],shapeDamage[-1]])
            damage1_layer = np.reshape(damage1_layer,[1,shapeDamage[0],shapeDamage[-1]])
            damage2_layer = np.reshape(damage2_layer,[1,shapeDamage[0],shapeDamage[-1]])
            damage3_layer = np.reshape(damage3_layer,[1,shapeDamage[0],shapeDamage[-1]])
            damage4_layer = np.reshape(damage4_layer,[1,shapeDamage[0],shapeDamage[-1]])
            damage5_layer = np.reshape(damage5_layer,[1,shapeDamage[0],shapeDamage[-1]])
            damage6_layer = np.reshape(damage6_layer,[1,shapeDamage[0],shapeDamage[-1]])
            damage7_layer = np.reshape(damage7_layer,[1,shapeDamage[0],shapeDamage[-1]])
            damage8_layer = np.reshape(damage8_layer,[1,shapeDamage[0],shapeDamage[-1]])



          # range0 = np.arange(0,int(nGauss_3/3),3)
          # Ranges for stresses and strains
          range0 = np.arange(0,int(nGauss_3),3)
          range1 = np.arange(1,int(nGauss_3),3)
          range2 = np.arange(2,int(nGauss_3),3)

          # Ranges for damage
          if use_damage:
            range0_damage = np.arange(0,int(nGauss_9),9)
            range1_damage = np.arange(1,int(nGauss_9),9)
            range2_damage = np.arange(2,int(nGauss_9),9)
            range3_damage = np.arange(3,int(nGauss_9),9)
            range4_damage = np.arange(4,int(nGauss_9),9)
            range5_damage = np.arange(5,int(nGauss_9),9)
            range6_damage = np.arange(6,int(nGauss_9),9)
            range7_damage = np.arange(7,int(nGauss_9),9)
            range8_damage = np.arange(8,int(nGauss_9),9)


          # Grids for stresses and strains
          grid0 = np.ix_([iLayer],range(numEle),range0)
          grid1 = np.ix_([iLayer],range(numEle),range1)
          grid2 = np.ix_([iLayer],range(numEle),range2)


          # Grids for damage
          if use_damage:
            grid0_damage = np.ix_([iLayer],range(numEle),range0_damage)
            grid1_damage = np.ix_([iLayer],range(numEle),range1_damage)
            grid2_damage = np.ix_([iLayer],range(numEle),range2_damage)
            grid3_damage = np.ix_([iLayer],range(numEle),range3_damage)
            grid4_damage = np.ix_([iLayer],range(numEle),range4_damage)
            grid5_damage = np.ix_([iLayer],range(numEle),range5_damage)
            grid6_damage = np.ix_([iLayer],range(numEle),range6_damage)
            grid7_damage = np.ix_([iLayer],range(numEle),range7_damage)
            grid8_damage = np.ix_([iLayer],range(numEle),range8_damage)


          # Assigning stresses and strains to the arrays
          sigma[grid0] = sigmaXX_layer
          sigma[grid1] = sigmaYY_layer
          sigma[grid2] = sigmaXY_layer
          strain[grid0] = strainXX_layer
          strain[grid1] = strainYY_layer
          strain[grid2] = strainXY_layer

          # Assigning damage variables to array
          if use_damage:
            damage[grid0_damage] = damage0_layer
            damage[grid1_damage] = damage1_layer
            damage[grid2_damage] = damage2_layer
            damage[grid3_damage] = damage3_layer
            damage[grid4_damage] = damage4_layer
            damage[grid5_damage] = damage5_layer
            damage[grid6_damage] = damage6_layer
            damage[grid7_damage] = damage7_layer
            damage[grid8_damage] = damage8_layer


        sigma_patches.append(sigma)
        strain_patches.append(strain)
        if use_damage:
            damage_patches.append(damage)


        prevNlayers = Nlayers


      # Until here I have the controlPts for each patch and the U vector for each patch (I hope)

    
      for (surf,i) in zip(surfList,range(len(surfList))):
        surf = surfList[i]

        shell={}

        Nlayers = Nlayers_patches[i]
        zLocs = zLocs_list[i]


        P = surf.degree_u
        Q = surf.degree_v

        noPtsX  = surf.ctrlpts_size_u
        noPtsY  = surf.ctrlpts_size_v


        uKnot = surf.knotvector_u
        vKnot = surf.knotvector_v

        # Recover info for the current patch from the lists
        controlPts = controlPts_patches[i]
        weights = [1]*len(controlPts)
        U = U_patches[i]
        sigma = sigma_patches[i]
        strain = strain_patches[i]
        if use_damage:
            damage = damage_patches[i]

        print('np.shape(sigma) = ', np.shape(sigma))


        shell['p'] = P
        shell['q'] = Q 
        shell['uKnot'] = uKnot 
        shell['vKnot'] = vKnot 
        shell['noPtsX'] = noPtsX 
        shell['noPtsY'] = noPtsY 
        # shell['weights'] = np.transpose(weights)
        shell['weights'] = weights
        shell['controlPts'] = controlPts
        shell['u'] = U 
        shell['stresses'] = sigma
        shell['strains'] = strain
        if use_damage:
            shell['damage'] = damage
        shell['nLayers'] = Nlayers
        shell['zLocs'] = zLocs


        matName = f'shell_{i}_{time}.mat'
        sio.savemat(matName,shell)

        print("Done saving .mat = ", matName)


exit()

# PostProcessing

import scipy.io as sio
import h5py

file = h5py.File('iga_cantilever.mpco', 'r')

patches = file['MODEL_STAGE[1]']['MODEL']['ELEMENTS']['1-UnknownMovableObject[0:0]'][:, 0] # These are the subdomains -> Iga Patches
elements = file['MODEL_STAGE[1]']['MODEL']['ELEMENTS']['207-UnknownMovableObject[0:0]'][:, 0] # These are the IGA_KLShell elements

nodes = file['MODEL_STAGE[1]']['MODEL']['ELEMENTS']['207-UnknownMovableObject[0:0]'][:, 1:]

print('nodes = \n', nodes)

numPatches = 1

nodesOnPatches = []
elementsOnPatches = []

continueSearch = True
while continueSearch:

  for i in range(1,len(nodes)):
    currentElement = elements[i-1]
    nextElement = elements[i]

    print("currentElement = ", currentElement)
    print("nextElement = ", nextElement)

    foundPatch = True

    if abs(nextElement - currentElement) > 1:
      continueSearch = True
      foundPatch = True
      numPatches += 1
      print('breaking in i = ', i)
      break_i = i
      nodesOnPatches.append(nodes[:i,:])
      elementsOnPatches.append(elements[:i])
      nodes = nodes[i:,:]
      elements = elements[i:]
      print('nodes = ', nodes)
      print('elements = ', elements)
      break
    else:
      foundPatch = False

  if not foundPatch:
    continueSearch = False
    nodesOnPatches.append(nodes)
    elementsOnPatches.append(elements)

  else:
    pass


# Getting nodeLimits in order to know which nodes are from which patch
nodeLimits = []
for nodes in nodesOnPatches:
  nodeLimits.append([nodes[0,0],nodes[-1,-1]])


nTimes = len(file['MODEL_STAGE[1]']['RESULTS']['ON_NODES']['DISPLACEMENT']['DATA'])

for time in range(nTimes):

  step = f'STEP_{time}'

  controlPts = file['MODEL_STAGE[1]']['MODEL']['NODES']['COORDINATES'][:, :] # All the control points / nodes in the model
  U = file['MODEL_STAGE[1]']['RESULTS']['ON_NODES']['DISPLACEMENT']['DATA'][step][:,:] # All the displacements in the model
  sigma_global = file['MODEL_STAGE[1]']['RESULTS']['ON_ELEMENTS']['section.fiber.stress']['207-UnknownMovableObject[0:0:0]']['DATA'][step][:,:] # Shape (numElements x Nlayers*noGps*3)
  strain_global = file['MODEL_STAGE[1]']['RESULTS']['ON_ELEMENTS']['section.fiber.strain']['207-UnknownMovableObject[0:0:0]']['DATA'][step][:,:] # Shape (numElements x Nlayers*noGps*3)
  zLocs = file['MODEL_STAGE[1]']['MODEL']['SECTION_ASSIGNMENTS']['SECTION_1[UnkownClassType]']['FIBER_DATA'][:,1] # Shape (numMaterial x 3 (algo, zLoc, thickness))

  controlPts_patches = []
  sigma_patches = []
  strain_patches =[]
  U_patches=[]

  for i in range(len(nodeLimits)):
    nodeStart, nodeEnd = nodeLimits[i]
    rangeNodes = np.arange(nodeStart,nodeEnd+1,1) - 1
    grid = np.ix_(rangeNodes,range(3))

    controlPts_here = controlPts[grid]
    U_here = U[grid]
    shapeU_here = np.shape(U_here)
    U_here = np.reshape(U_here,shapeU_here[0]*shapeU_here[1])


    elementsThisPatch = elementsOnPatches[i]

    controlPts_patches.append(controlPts_here)
    U_patches.append(U_here)

    print('elementsThisPatch = ', elementsThisPatch)

    # Have to do a new grid
    if i == 0:
      # I'm in the first patch
      numEleThisPatch = len(elementsThisPatch)
      id_elements = np.arange(0,numEleThisPatch,1)
    else:
      elementsPrevPatch = elementsOnPatches[i-1]
      numEleThisPatch = len(elementsThisPatch)
      numElePrevPatch = len(elementsPrevPatch)
      id_elements = np.arange(numElePrevPatch, numElePrevPatch+numEleThisPatch,1)

    print('id_elements = ', id_elements)

    nData = len(sigma_global[0])

    grid = np.ix_(id_elements,range(nData))

    # continue

    # Getting stresses and strains from elements in each layer
    sigma = sigma_global[grid]
    strain = strain_global[grid]

    print('sigma = ', sigma)


    shapeSigma=np.shape(sigma) # numEle x 3*nGauss
    shapeStrain=np.shape(strain)

    numEle = shapeSigma[0] 

    # print('numEle = ', numEle)

    range0 = np.arange(0,shapeSigma[1],3)
    range1 = range0+1
    range2 = range0+2

    print('range0 = ', range0)

    # Getting quantities in each direction for all layers
    sigmaXX = sigma[:,np.ix_(range0)]
    sigmaYY = sigma[:,np.ix_(range1)]
    sigmaXY = sigma[:,np.ix_(range2)]
    strainXX = strain[:,np.ix_(range0)]
    strainYY = strain[:,np.ix_(range1)]
    strainXY = strain[:,np.ix_(range2)]


    print('sigmaXX = ', sigmaXX)

    nGauss_3 = int(shapeSigma[1]/Nlayers)

    print('nGauss_3 = ', nGauss_3)

    sigma = np.zeros([Nlayers,numEle,nGauss_3])
    strain = np.zeros([Nlayers,numEle,nGauss_3])



    for iLayer in range(Nlayers):

      rangeLayer = np.arange(iLayer,int(nGauss_3/3)*Nlayers,Nlayers)
      # if Nlayers == 1:
      #   rangeLayer = np.arange(iLayer,int(nGauss_3/3),Nlayers)
      # else:
      #   rangeLayer = np.arange(iLayer,nGauss_3,Nlayers)

      print('rangeLayer = ', rangeLayer)

      sigmaXX_layer = sigmaXX[:,:,np.ix_(rangeLayer)]
      sigmaYY_layer = sigmaYY[:,:,np.ix_(rangeLayer)]
      sigmaXY_layer = sigmaXY[:,:,np.ix_(rangeLayer)]

      strainXX_layer = strainXX[:,:,np.ix_(rangeLayer)]
      strainYY_layer = strainYY[:,:,np.ix_(rangeLayer)]
      strainXY_layer = strainXY[:,:,np.ix_(rangeLayer)]

      shapeSigma = np.shape(sigmaXX_layer)
      shapeStrain = np.shape(strainXX_layer)

      # print('shapeSigma = ', shapeSigma)

      sigmaXX_layer = np.reshape(sigmaXX_layer,[1,shapeSigma[0],shapeSigma[-1]])
      sigmaYY_layer = np.reshape(sigmaYY_layer,[1,shapeSigma[0],shapeSigma[-1]])
      sigmaXY_layer = np.reshape(sigmaXY_layer,[1,shapeSigma[0],shapeSigma[-1]])

      strainXX_layer = np.reshape(strainXX_layer,[1,shapeStrain[0],shapeStrain[-1]])
      strainYY_layer = np.reshape(strainYY_layer,[1,shapeStrain[0],shapeStrain[-1]])
      strainXY_layer = np.reshape(strainXY_layer,[1,shapeStrain[0],shapeStrain[-1]])


      # range0 = np.arange(0,int(nGauss_3/3),3)
      range0 = np.arange(0,int(nGauss_3),3)
      range1 = np.arange(1,int(nGauss_3),3)
      range2 = np.arange(2,int(nGauss_3),3)

      grid0 = np.ix_([iLayer],range(numEle),range0)
      grid1 = np.ix_([iLayer],range(numEle),range1)
      grid2 = np.ix_([iLayer],range(numEle),range2)


      sigma[grid0] = sigmaXX_layer
      sigma[grid1] = sigmaYY_layer
      sigma[grid2] = sigmaXY_layer


      strain[grid0] = strainXX_layer
      strain[grid1] = strainYY_layer
      strain[grid2] = strainXY_layer


      sigma_patches.append(sigma)
      strain_patches.append(strain)


  # Until here I have the controlPts for each patch and the U vector for each patch (I hope)

  surfList = [surfL, surfR]

  for (surf,i) in zip(surfList,range(len(surfList))):

    shell={}


    P =surf.degree_u
    Q =surf.degree_v

    noPtsX  = surf.ctrlpts_size_u
    noPtsY  = surf.ctrlpts_size_v


    uKnot = surf.knotvector_u
    vKnot = surf.knotvector_v

    # Recover info for the current patch from the lists
    controlPts = controlPts_patches[i]
    weights = [1]*len(controlPts)
    U = U_patches[i]
    sigma = sigma_patches[i]
    strain = strain_patches[i]


    shell['p'] = P
    shell['q'] = Q 
    shell['uKnot'] = uKnot 
    shell['vKnot'] = vKnot 
    shell['noPtsX'] = noPtsX 
    shell['noPtsY'] = noPtsY 
    # shell['weights'] = np.transpose(weights)
    shell['weights'] = weights
    shell['controlPts'] = controlPts
    shell['u'] = U 
    shell['stresses'] = sigma
    shell['strains'] = strain
    shell['nLayers'] = Nlayers
    shell['zLocs'] = zLocs

    matName = f'shell_{i}_{time}.mat'
    sio.savemat(matName,shell)

    print("Done saving .mat")