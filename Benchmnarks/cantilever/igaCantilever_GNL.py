#  IGA CANTILEVER PLATE with bending strip


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
  [0.00000  , 1.0 , 0.0 , 1.0] ,
  [0.31250  , 0.0 , 0.0 , 1.0] ,
  [0.31250  , 1.0 , 0.0 , 1.0] ,
  [0.93750  , 0.0 , 0.0 , 1.0] ,
  [0.93750  , 1.0 , 0.0 , 1.0] ,
  [1.87500  , 0.0 , 0.0 , 1.0] ,
  [1.87500  , 1.0 , 0.0 , 1.0] ,
  [3.12500  , 0.0 , 0.0 , 1.0] ,
  [3.12500  , 1.0 , 0.0 , 1.0] ,
  [4.37500  , 0.0 , 0.0 , 1.0] ,
  [4.37500  , 1.0 , 0.0 , 1.0] ,
  [5.62500  , 0.0 , 0.0 , 1.0] ,
  [5.62500  , 1.0 , 0.0 , 1.0] ,
  [6.87500  , 0.0 , 0.0 , 1.0] ,
  [6.87500  , 1.0 , 0.0 , 1.0] ,
  [8.12500  , 0.0 , 0.0 , 1.0] ,
  [8.12500  , 1.0 , 0.0 , 1.0] ,
  [9.06250  , 0.0 , 0.0 , 1.0] ,
  [9.06250  , 1.0 , 0.0 , 1.0] ,
  [9.68750  , 0.0 , 0.0 , 1.0] ,
  [9.68750  , 1.0 , 0.0 , 1.0] ,
  [10.00000 , 0.0 , 0.0 , 1.0] ,
  [10.00000 , 1.0 , 0.0 , 1.0]
])

weights =  controlPts[:,-1]


# Create a NURBS surface instance
surf = NURBS.Surface()

# Set surface degrees
surf.degree_u = 4
surf.degree_v = 1

# Setting control points for surface
surf.set_ctrlpts(controlPts.tolist(), 12, 2)

# Set knot vectors
surf.knotvector_u = knotvector.generate(surf.degree_u, surf.ctrlpts_size_u)
surf.knotvector_v = knotvector.generate(surf.degree_v, surf.ctrlpts_size_v)


# nDMaterial ElasticIsotropic $nDtag_elastic $elasticidad_probeta
# $poisson_probeta
E1 = 1.2e6  # Young's modulus N/m^2
E2 = E1
nu12 = 0.0  # Poisson's ratio
nu21 = 0.0  # Poisson's ratio
G12 = E1/(2*(1+nu21)) # Shear modulus
rho = 8.0e3  # *9.807 # kg/m^3
t = 100. * mm


tagPlaneStress1 = 1
ops.nDMaterial("ElasticOrthotropicPlaneStress", tagPlaneStress1, E1, E2, nu12, nu21, G12, rho)



deg2rad = pi / 180

matTags = [1,1,1]
thickness = [t/3,t/3,t/3]
θ = [0 * deg2rad,0 * deg2rad, 0 * deg2rad]

gFact = [0.0, 0.0, 0.0 * 9.807]


Nlayers = len(θ)

patchTag=1
nodeStartTag=1


# Flipping control point to u,v format
controlPts = surf.ctrlpts2d[:]
controlPts = np.array(compatibility.flip_ctrlpts2d(controlPts))

for dim in controlPts: # Unweighting control pts
    for point in dim:
        point[0:3]/=point[3]

# Creating a Patch in OpenSees
ops.IGA("Patch", patchTag, nodeStartTag, surf.degree_u, surf.degree_v, surf.ctrlpts_size_u, surf.ctrlpts_size_v,
        "-type", "KLShell",
        # "-nonLinearGeometry", 0,
        "-planeStressMatTags", *matTags,
        "-gFact", *gFact,
        "-theta", *θ,
        "-thickness", *thickness,
        "-uKnot", *surf.knotvector_u, "-vKnot", *surf.knotvector_v, "-controlPts", *controlPts.flatten())




# Creating constraints
for n in ops.getNodeTags():
    n = int(n)
    if n in [1, 2, 13, 14]:
        ops.fix(n, 1, 1, 1)


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
for n in [12, 24]:
    # # Vertical load
    ops.load(n, 0, 0, Pz / 2.0)

I = Lb * (sum(thickness)**3) / 12.0

print("Finished loading nodes")


print("Starting analysis")

# Create test
ops.test("NormDispIncr", 1.0e-8, 500, 1)


# create SOE
ops.system("UmfPack")

# create DOF number
ops.numberer("Plain")

# create constraint handler
ops.constraints("Plain")
# ops.constraints("Penalty",1,1)


# ops.algorithm("Linear")
# ops.algorithm("Newton")
# ops.algorithm("SecantNewton")
# ops.algorithm("NewtonLineSearch")
# ops.algorithm("ModifiedNewton")
ops.algorithm("KrylovNewton")
# ops.algorithm("BFGS")
# ops.algorithm("Broyden")

# create integrator

# nSteps = 50
# ops.integrator("LoadControl", 1.0 / nSteps)

delta = 0.67/2
defMax = 6.7
nSteps = abs(int(defMax / delta))
# ops.integrator("LoadControl", 1.0 / nSteps)
ops.integrator("DisplacementControl", 12, 3, delta)


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
    data[j + 1, 0] = abs(ops.nodeDisp(12, 3))
    data[j + 1, 1] = abs(ops.nodeDisp(12, 1))
    data[j + 1, 2] = abs(ops.getLoadFactor(1) * Pz)
    print(data[j + 1, 2])





# Adding deformation to control points
fdef = 1e0
nodes=ops.getNodeTags()
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



# Visualize surface
surfVisualize(surf, hold=False)


result = np.array(ops.nodeDisp(24, 3))
print("result: ", result)

F=Pz*np.arange(0,1.05,0.05)
U=[0,0.026,0.103,0.224,0.381,0.563,0.763,0.971,1.184,1.396,1.604,1.807,2.002,2.19,2.37,2.541,2.705,2.861,3.01,3.151,3.286]
W=[0,0.663,1.309,1.922,2.493,3.015,3.488,3.912,4.292,4.631,4.933,5.202,5.444,5.66,5.855,6.031,6.19,6.335,6.467,6.588,6.698]
plt.plot(data[:, 0], data[:, 2], 'or', mfc="None", label="W disp (IGA)") #w
plt.plot(W, F, '-r',label="W disp (reference)") #w

plt.plot(data[:, 1], data[:, 2], 'ob', mfc="None", label="U disp (GNL)") #u
plt.plot(U, F, '-b',label="U disp (reference)") #u

plt.xlabel('Horizontal Displacement')
plt.ylabel('Vertical Load')
plt.legend()
plt.show()

print("Done")


print("Done!")

ops.remove("recorders")


import scipy.io as sio
import h5py

file = h5py.File('iga_cantilever.mpco', 'r')

patches = file['MODEL_STAGE[1]']['MODEL']['ELEMENTS']['1-UnknownMovableObject[0:0]'][:, 0]
elements = file['MODEL_STAGE[1]']['MODEL']['ELEMENTS']['207-UnknownMovableObject[0:0]'][:, 0]

print('patches = ', patches)
print('elements = ', elements)
# exit()

controlPts = file['MODEL_STAGE[1]']['MODEL']['NODES']['COORDINATES'][:, :]
# print(controlPts)
# controlPts = np.array(compatibility.flip_ctrlpts(controlPts,noPtsY,noPtsX))
# print(controlPts)
# exit()

U = file['MODEL_STAGE[1]']['RESULTS']['ON_NODES']['DISPLACEMENT']['DATA']['STEP_19'][:,:]
shapeU = np.shape(U)
U = np.reshape(U,shapeU[0]*shapeU[1])

sigma = file['MODEL_STAGE[1]']['RESULTS']['ON_ELEMENTS']['section.fiber.stress']['207-UnknownMovableObject[0:0:0]']['DATA']['STEP_19'][:,:] # Shape (numElements x nLayers*noGps*3)
strain = file['MODEL_STAGE[1]']['RESULTS']['ON_ELEMENTS']['section.fiber.strain']['207-UnknownMovableObject[0:0:0]']['DATA']['STEP_19'][:,:] # Shape (numElements x nLayers*noGps*3)



shapeSigma=np.shape(sigma) # numEle x 3*nGauss
shapeStrain=np.shape(strain)

numEle = shapeSigma[0] 

# print('numEle = ', numEle)

range0 = np.arange(0,shapeSigma[1],3)
range1 = range0+1
range2 = range0+2

# print('range0 = ', range0)

sigmaXX = sigma[:,np.ix_(range0)]
sigmaYY = sigma[:,np.ix_(range1)]
sigmaXY = sigma[:,np.ix_(range2)]

strainXX = strain[:,np.ix_(range0)]
strainYY = strain[:,np.ix_(range1)]
strainXY = strain[:,np.ix_(range2)]

nGauss_3 = int(shapeSigma[1]/Nlayers)

sigma = np.zeros([Nlayers,numEle,nGauss_3])
strain = np.zeros([Nlayers,numEle,nGauss_3])



for iLayer in range(Nlayers):

  # rangeLayer = np.arange(iLayer,len(sigmaXX[0,0]),2)
  rangeLayer = np.arange(iLayer,nGauss_3,Nlayers)

  # print('rangeLayer = ', rangeLayer)

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



zLocs = file['MODEL_STAGE[1]']['MODEL']['SECTION_ASSIGNMENTS']['SECTION_1[UnkownClassType]']['FIBER_DATA'][:,1] # Shape (numMaterial x 3 (algo, zLoc, thickness))

shell={}


P =surf.degree_u
Q =surf.degree_v

noPtsX  = surf.ctrlpts_size_u
noPtsY  = surf.ctrlpts_size_v


uKnot = surf.knotvector_u
vKnot = surf.knotvector_v


shell['p'] = P
shell['q'] = Q 
shell['uKnot'] = uKnot 
shell['vKnot'] = vKnot 
shell['noPtsX'] = noPtsX 
shell['noPtsY'] = noPtsY 
# shell['weights'] = np.transpose(weights)
shell['weights'] = weights
shell['controlPts'] = controlPts[:,0:3]
shell['u'] = U 
shell['stresses'] = sigma
shell['strains'] = strain
shell['nLayers'] = Nlayers
shell['zLocs'] = zLocs


sio.savemat('shell.mat',shell)

print("Done saving .mat")