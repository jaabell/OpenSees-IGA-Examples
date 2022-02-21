
#  IGA Scordelis Lo Roof
import numpy as np
import opensees as ops
from math import *

# Geomgl utilities for visualization and surface manipulation
from geomdl import NURBS, compatibility, operations, knotvector, multi
from surfVisualize import *


from oct2py import octave

octave.addpath(octave.genpath('octaveFiles/'))

mm = 1.0 / 1000.  # m

ops.wipe()
ops.model('basic', '-ndm', 3, '-ndf', 3)

refineLevel = 5
refU=3
refV=3
octave.makeScordelis(refineLevel,refU,refV)
# out = octave.makeCantileverShell()

import scipy.io as sio
shell = sio.loadmat('shell.mat')
import numpy

P = int(shell['p'])
Q = int(shell['q'])
uKnot = shell['uKnot'][0]
vKnot = shell['vKnot'][0]
noPtsX = int(shell['noPtsX'])
noPtsY = int(shell['noPtsY'])
weights = shell['weights']
controlPts = shell['controlPts']
for point in controlPts: #weighting controlPts
    point[0:3]*=point[3]




print('noPtsX = ', noPtsX)
print('noPtsY = ', noPtsY)
print('uKnot = ', uKnot)
print('vKnot = ', vKnot)



# exit()




patchTag = 1
P = P
Q = Q


# Create a NURBS surface instance
surf = NURBS.Surface()

# Set surface degrees
surf.degree_u = P
surf.degree_v = Q

# Setting control points for surface
surf.set_ctrlpts(controlPts.tolist(), noPtsY, noPtsX)
controlPts = surf.ctrlpts2d[:]
controlPts = compatibility.flip_ctrlpts2d(controlPts)
controlPts=numpy.reshape(controlPts,[noPtsX*noPtsY,4])
surf.set_ctrlpts(controlPts.tolist(), noPtsX, noPtsY)

# Set knot vectors
surf.knotvector_u = knotvector.generate(surf.degree_u, surf.ctrlpts_size_u)
surf.knotvector_v = knotvector.generate(surf.degree_v, surf.ctrlpts_size_v)

# Visualize surface
surfVisualize(surf, hold=True)




# nDMaterial ElasticIsotropic $nDtag_elastic $elasticidad_probeta
# $poisson_probeta
E1 = 4.32e8  # Young's modulus N/m^2
E2 = E1
nu = 0.0  # Poisson's ratio
g = -9.8066  # gravity, m/s^2
t = 0.25
rho = 90/abs(g)/t  # kg/m^3


tagNDmat1 = 1
ops.nDMaterial("ElasticIsotropic", tagNDmat1, E1, nu, rho)

tagNDmat2 = 2
ops.nDMaterial("ElasticIsotropic", tagNDmat2, E2, nu, rho)


# nDMaterial PlateFiber $nDtag_platefiber $nDtag_elastic
tagPlaneStress1 = 3
ops.nDMaterial("PlaneStress", tagPlaneStress1, tagNDmat1)

tagPlaneStress2 = 4
ops.nDMaterial("PlaneStress", tagPlaneStress2, tagNDmat2)

deg2rad = pi / 180

gFact = [0.0, 0.0, 0.0]

matTags = [3]
thickness = [t]
θ = [0 * deg2rad]


Nlayers = len(θ)

controlPts = surf.ctrlpts2d[:]
controlPts = np.array(compatibility.flip_ctrlpts2d(controlPts))

for point in controlPts:
    for i in range(surf.ctrlpts_size_v):
        for k in range(3):
            point[i][k] /= point[i][3]

nodeStartTag=1

ops.IGA("SurfacePatch", patchTag, nodeStartTag, surf.degree_u, surf.degree_v, surf.ctrlpts_size_u, surf.ctrlpts_size_v,
        "-type", "KLShell",
        "-nonLinearGeometry", 0,
        "-planeStressMatTags", *matTags,
        "-gFact", *gFact,
        "-theta", *θ,
        "-thickness", *thickness,
        "-uKnot", *surf.knotvector_u, "-vKnot", *surf.knotvector_v, "-controlPts", *controlPts.flatten())


print("\n\n\nPRINTING DOMAIN-----------------------")
ops.printModel()
print("\n\n\nDONE PRINTING DOMAIN-----------------------")


#  Dirichlet BCs (symmetry conditions)
#    z
#    |
#  (D)------- (C)
#    |      |
#    |    L |
#    |  R   |
#  (A)------- (B) --->x
#
#  rigid diaphram: AB: u_x, u_y = 0
#  Symmetry conditions on CD and AD
#  BC is free

nPoints = surf.ctrlpts_size_v * surf.ctrlpts_size_u
nodesOnCB = np.arange(1, nPoints, surf.ctrlpts_size_v)
nodesOnAD = nodesOnCB + surf.ctrlpts_size_v - 1
nodesOnAB = np.arange(1, surf.ctrlpts_size_u + 1)
nodesOnDC = np.arange(nPoints - surf.ctrlpts_size_u + 1, nPoints + 1)

nodesNextToAB = nodesOnAB + surf.ctrlpts_size_u
nodesNextToAD = nodesOnAD - 1

recordNode = 1 # Measure deformation here


for n in ops.getNodeTags():
    n = int(n)
    if n in nodesOnAB:
        ops.fix(n, 0, 0, 1)
    elif n in nodesOnAD:
        ops.fix(n, 1, 0, 0)
    elif n in nodesOnDC:
        ops.fix(n, 1, 1, 0)

# Symmetry conditions
for i in range(len(nodesOnDC)):
    # ops.equalDOF(int(nodesOnDC[i]), int(nodesNextToDC[i]), 1, 2, 3)
    # ops.equalDOF(int(nodesOnAD[i]), int(nodesNextToAD[i]), 1, 2, 3)

    ops.equalDOF(int(nodesOnAB[i]), int(nodesNextToAB[i]), 2)
    ops.equalDOF(int(nodesOnAD[i]), int(nodesNextToAD[i]), 2)

    # ops.equalDOF(int(nodesOnAB[i]), int(nodesNextToAB[i]), 1,2,3)
    # ops.equalDOF(int(nodesOnAD[i]), int(nodesNextToAD[i]), 1,2,3)


# ------------------------------
# Start of analysis generation
# ------------------------------

# create TimeSeries
ops.timeSeries("Linear", 1)

# create a plain load pattern
ops.pattern("Plain", 1, 1)


print("Loading nodes")
weight = [0.0, g, 0.0]
ops.eleLoad("-ele", 1, "-type", "-SelfWeight", *weight)
print("Finished loading nodes")


# create SOE
# ops.system("FullGeneral")
ops.system("UmfPack")

# ops.system("BandSPD")

# create DOF number
ops.numberer("Plain")


print("Starting analysis")

# create constraint handler
ops.constraints("Plain")
# ops.constraints("Penalty", 1e12*E1,1e12*E1)

# Create test
ops.test("NormDispIncr", 1.0e-9, 50, 1)
# ops.test("NormUnbalance",1e-8,10)

# create algorithm

ops.algorithm("Linear")
# ops.algorithm("Newton")
# ops.algorithm("ModifiedNewton")
# ops.algorithm("KrylovNewton")

# create integrator
# nSteps=10
# ops.integrator("LoadControl", 1.0/nSteps)
ops.integrator("LoadControl", 1.0)


# create analysis object
ops.analysis("Static")

ops.analyze(1)

print("Finished analysis")


# Adding deformation to controlPts
controlPts = surf.ctrlpts2d[:]
controlPts = compatibility.flip_ctrlpts2d(controlPts)  # Flipping to u,v

fDefX = 1e1
fDefY = 1e1
fDefZ = 1e1
i = 1
for dim in controlPts:
    for point in dim:
        weight=point[3]

        point[0]/=weight
        point[1]/=weight
        point[2]/=weight

        point[0] = (ops.nodeCoord(i,1)+fDefX * ops.nodeDisp(i,1)) * weight  # Times the weight
        point[1] = (ops.nodeCoord(i,2)+fDefY * ops.nodeDisp(i,2)) * weight  # Times the weight
        point[2] = (ops.nodeCoord(i,3)+fDefZ * ops.nodeDisp(i,3)) * weight  # Times the weight

        i += 1

# Setting control points for surface
controlPts = compatibility.flip_ctrlpts2d(controlPts)
controlPts = (np.array(controlPts).reshape(
    nPoints, 4))
surf.set_ctrlpts(controlPts.tolist(), surf.ctrlpts_size_u, surf.ctrlpts_size_v)



# Symmetry
bbox = np.array(surf.bbox)
midPoint = (bbox[1] - bbox[0]) / 2


surfSym = operations.rotate(surf, 0, axis=1, inplace=False)
ctrlpts_surfSym = surfSym.ctrlpts2d[:]
for dim in ctrlpts_surfSym:
    for point in dim:
        point[0] *= -1  # Symmetry along x axis

ctrlpts_surfSym.reverse()  # Flipping
ctrlpts_surfSym = (np.array(ctrlpts_surfSym).reshape(surfSym.ctrlpts_size_u * surfSym.ctrlpts_size_v, 4))
surfSym.set_ctrlpts(ctrlpts_surfSym.tolist(), surfSym.ctrlpts_size_u, surfSym.ctrlpts_size_v)



surfSym_01=operations.rotate(surf, 0, axis=1, inplace=False)
ctrlpts_surfSym_01 = surfSym_01.ctrlpts2d[:]
for dim in ctrlpts_surfSym_01:
    for point in dim:
        point[2] *= -1  # Symmetry along x axis

ctrlpts_surfSym_01.reverse()  # Flipping
ctrlpts_surfSym_01 = (np.array(ctrlpts_surfSym_01).reshape(surfSym_01.ctrlpts_size_u * surfSym_01.ctrlpts_size_v, 4))
surfSym_01.set_ctrlpts(ctrlpts_surfSym_01.tolist(), surfSym_01.ctrlpts_size_u, surfSym_01.ctrlpts_size_v)

surfSym_02=operations.rotate(surfSym, 0, axis=1, inplace=False)
ctrlpts_surfSym_02 = surfSym_02.ctrlpts2d[:]
for dim in ctrlpts_surfSym_02:
    for point in dim:
        point[2] *= -1  # Symmetry along x axis

ctrlpts_surfSym_02.reverse()  # Flipping
ctrlpts_surfSym_02 = (np.array(ctrlpts_surfSym_02).reshape(surfSym_02.ctrlpts_size_u * surfSym_02.ctrlpts_size_v, 4))
surfSym_02.set_ctrlpts(ctrlpts_surfSym_02.tolist(), surfSym_02.ctrlpts_size_u, surfSym_02.ctrlpts_size_v)

# Creating container for multipatches

surfList = [surf, surfSym, surfSym_01, surfSym_02]


container = multi.SurfaceContainer(surfList)

# Visualize surface

container.sample_size = 30
for surf in container:
    surf.evaluate()

# Visualization configuration
container.vis = VisVTK.VisSurface(ctrlpts=False, legend=False, line_width=1, trim_size=20)
# container.vis.ctrlpts_offset=0.1

# Render the surface
evalcolor = ["red", "green", "green", "green"]
cpcolor=["red","black", "black", "black"]
container.render(evalcolor=evalcolor, cpcolor=cpcolor)


# # Visualize surface
# surfVisualize(surf, hold=True)

dispC=ops.nodeDisp(recordNode)
print("dispC: ", dispC)


print("Done")


