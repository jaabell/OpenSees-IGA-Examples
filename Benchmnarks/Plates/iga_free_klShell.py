# free thin plate verification using Isogeometric Kirchhoff Love Shell
#  python using geomdl
# Felipe Elgueta
# Universidad de los Andes, Chile

import os
import sys
import numpy as np
# from geomdl import BSpline
from geomdl import NURBS
from geomdl import operations
from geomdl import exchange
from geomdl import utilities
from geomdl import helpers
from geomdl import compatibility, knotvector
from geomdl.visualization import VisVTK as vis

from scipy.sparse import csr_matrix, lil_matrix, linalg
from scipy.sparse.linalg import eigsh


from surfVisualize import *

import opensees as ops


np.set_printoptions(formatter={'float': lambda x: "{0:0.6f}".format(x)})

mm=1.0/1000.0

La = 341.4*mm  # m
Lb = 318.8*mm  # m
t = 3.53*mm  # m

controlPts = np.array([
    [-La / 2    , -Lb / 2    , 0 , 1] ,
    [-La / 2    , 0          , 0 , 1] ,
    [-La / 2    , Lb / 2     , 0 , 1] ,
    [0          , -Lb / 2    , 0 , 1] ,
    [0          , 0          , 0 , 1] ,
    [0          , Lb / 2     , 0 , 1] ,
    [La / 2     , -Lb / 2    , 0 , 1] ,
    [La / 2     , 0          , 0 , 1] ,
    [La / 2     , Lb / 2     , 0 , 1]
])
p = 2
q = 2


# Create a BSpline surface instance
surf = NURBS.Surface()


# Set degrees
surf.degree_u = p 
surf.degree_v = q 

# Set control points
surf.set_ctrlpts(controlPts.tolist(), 3, 3)

# Set knot vectors
surf.knotvector_u = knotvector.generate(surf.degree_u, surf.ctrlpts_size_u)
surf.knotvector_v = knotvector.generate(surf.degree_v, surf.ctrlpts_size_v)


# Refine knot vectors and update geometry (refine both directions)
# u,v
refU = 3
refV = 3
operations.refine_knotvector(surf, [refU, refV])

# Get new control points and weights after refining
weights = surf.weights
controlPts = surf.ctrlpts2d[:]

# Update geometry after refining
noPtsX = surf.ctrlpts_size_u
noPtsY = surf.ctrlpts_size_v


# Set degrees
if refU == 0:
    p = 2
if refV == 0:
    q = 2

p=5
q=5
surf.degree_u = p
surf.degree_v = q

surf.knotvector_u = knotvector.generate(surf.degree_u, surf.ctrlpts_size_u)
surf.knotvector_v = knotvector.generate(surf.degree_v, surf.ctrlpts_size_v)

noCtrPts = len(controlPts)

# Visualize surface
surfVisualize(surf, hold=True)


# Material parameters
GPa=1e9
E1 = 152.7*GPa # Young's modulus N/m^2
E2 = 8.832*GPa
nu12 = 0.297  # Poisson's ratio
nu21 = nu12*E2/E1  # Poisson's ratio
G12=5.274*GPa
rho = 1560  # kg/m^3



tagPlaneStress1 = 1
ops.nDMaterial("ElasticOrthotropicPlaneStress", tagPlaneStress1, E1, E2, nu12, nu21, G12, rho)



#  Dirichlet BCs (symmetry conditions)
#    z
#    |
#  (D)------- (C)
#    |      |
#    |    L |
#    |  R   |
#  (A)------- (B) --->x
#
#  AB: free
#  Symmetry conditions on BC,CD and AD
#

nPoints = surf.ctrlpts_size_u * surf.ctrlpts_size_v

nodesOnDC = np.arange(1, surf.ctrlpts_size_u + 1, 1)
nodesNextToDC = nodesOnDC + surf.ctrlpts_size_u

nodesOnCB = np.arange(1, nPoints, surf.ctrlpts_size_v)
nodesNextToCB = nodesOnCB + 1


nodesOnAD = np.arange(surf.ctrlpts_size_u, nPoints + 1, surf.ctrlpts_size_v)
nodesNextToAD = nodesOnAD - 1


deg2rad = np.pi / 180

# matTags = [3, 4, 3, 4, 3]
# thickness = [10. * mm, 10. * mm, 10. * mm, 10. * mm, 10. * mm]
# θ = [0 * deg2rad, 45 * deg2rad, 90 * deg2rad, -45 * deg2rad, 0 * deg2rad]

θ = [-15,15,15,-15,15,-15,-15,15,-15,15,15,-15]
θ = (np.array(θ)*2*deg2rad).tolist()
Nlayers = len(θ)
thickness = [t/Nlayers]*Nlayers
matTags = [1]*len(θ)

print("thickness: ", thickness)


gFact = [0.0, 0.0, 0.0]



shellType = 'KLShell'
patchTag = 1
nodeStartTag = 1

nodesMap = []
patchTags=[]

patchTags.append(patchTag)
    

# Flipping control point to u,v format
controlPts = surf.ctrlpts2d[:]
controlPts = np.array(compatibility.flip_ctrlpts2d(controlPts))

for dim in controlPts:  # Unweighting control pts
    for point in dim:
        point[0:3] /= point[3]

# Creating a Patch in OpenSees
ops.IGA("Patch", patchTag, nodeStartTag, surf.degree_u, surf.degree_v, surf.ctrlpts_size_u, surf.ctrlpts_size_v,
        "-type", shellType,
        "-nonLinearGeometry", 0,
        "-planeStressMatTags", *matTags,
        "-gFact", *gFact,
        "-theta", *θ,
        "-thickness", *thickness,
        "-uKnot", *surf.knotvector_u, "-vKnot", *surf.knotvector_v, "-controlPts", *controlPts.flatten())

# Get the nodes on current patch
nodesMap.append(np.arange(nodeStartTag, ops.getNodeTags()[-1] + 1).tolist())

# Update patchTag, nodeStartTag and materialTag
lastElTag = ops.getEleTags()[-1]
lastNodeTag = ops.getNodeTags()[-1]
matTags[0] += 1
patchTag = lastElTag + 1
nodeStartTag = lastNodeTag + 1


fixedNodes=np.concatenate([nodesOnCB,nodesNextToCB])
# nextToFixedNodes=nodesNextToAD

for n in fixedNodes:
    n=int(n)
    ops.fix(n,1,1,1)


print("\n\n\nPRINTING DOMAIN-----------------------")
ops.printModel()
print("\n\n\nDONE PRINTING DOMAIN-----------------------")


nodes = ops.getNodeTags()

Nnodes = len(nodes)
Neigenvalues = 7  # arpack can only compute N-1 eigvals


w2s = ops.eigen(Neigenvalues)
# w2s = ops.eigen('-solver','-fullGenLapack',Neigenvalues)
# w2s = ops.eigen('-standard','-symmBandLapack',Neigenvalues)


order = np.argsort(w2s)
w2s = np.array(w2s, dtype=np.float64)[order]
w=np.sqrt(w2s)
print("w: ", w)
# exit()
# exit()


for i, w2 in enumerate(w2s):
    w = np.sqrt(abs(w2))
    f = w / 2 / np.pi
    T = 1 / f
    print(f"{i} {w2} {w} {f} {T} ")


phi = np.zeros((3 * Nnodes, Neigenvalues))

for i, n in enumerate(nodes):
    for j in range(Neigenvalues):
        phi[3 * i:3 * i + 3, j] = ops.nodeEigenvector(n, j + 1)

print(f"ϕ = {phi}")

w = np.sqrt(abs(w2s))
print("w/2/pi: ", w/2/np.pi)
for j in range(Neigenvalues):
    w = np.sqrt(abs(w2s[j]))
    f = w / 2 / np.pi
    print("=================================")
    print(f"Eigenvalue {j}")
    print(f"f = {f} Hz")
    print("=================================")


    # Adding deformation to control points
    fdef = 0.5e-1

    nodes = nodesMap[0]
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
    # surfVisualize(surf, hold=False)

    
    nodes = nodesMap[0]
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