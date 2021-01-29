#  IGA CANTILEVER PLATE. ONE ELEMENT MESH, LINEAR CONVERGENCE OBTAINED


import numpy as np
import opensees as ops
from math import *

# Geomgl utilities for visualization and surface manipulation
from geomdl import NURBS, compatibility, operations
from surfVisualize import *


def getCtrlPtsAndWeights(surf):
    # Flip control points to [u][v] (they are given by the getters in [v][u])
    noPtsX = surf.ctrlpts_size_u
    noPtsY = surf.ctrlpts_size_v
    weights = surf.weights
    controlPts = surf.ctrlpts2d[:]
    controlPts = compatibility.flip_ctrlpts2d(controlPts)
    controlPts = (np.array(controlPts).reshape(noPtsX * noPtsY, 4)).tolist()
    noCtrPts = len(controlPts)

    # Get separate array of control points in u,v and weights

    for i in range(len(controlPts)):
        pt = controlPts[i][:]
        wt = pt[-1]
        pt[0:3] = np.array(pt[0:3]) / wt
        controlPts[i] = pt[0:3]
        weights[i] = wt

    return controlPts, weights


def generateKnotVector(deg, nPts):

    import numpy as np
    knotVector = np.zeros(nPts + deg + 1)
    nMiddle = len(knotVector) - 2 * (deg + 1)
    step = 1.0 / (nMiddle + 1)

    ini = np.zeros(deg + 1)
    if step == 0.5:
        middle = np.array([step])
    else:
        middle = np.arange(0 + step, 1 - step, step)
        middle = np.linspace(0 + step, 1 - step, nPts +
                             deg + 1 - 2 * (deg + 1))
    fin = ini + 1

    knotVector = np.copy(ini)
    knotVector = np.append(knotVector, middle)
    knotVector = np.append(knotVector, fin)

    return knotVector


La = 10.0  	#
Lb = 1.0  	#
mm = 1.0 / 1000.  # m

ops.wipe()
ops.model('basic', '-ndm', 3, '-ndf', 3)

# These are given in v,u
controlPts = np.array([
    [0, 0, 0, 1],
    [La * 1 / 3, 0, 0, 1],
    [La * 2 / 3, 0, 0, 1],
    [La, 0, 0, 1],
    [0, Lb, 0, 1],
    [La * 1 / 3, Lb, 0, 1],
    [La * 2 / 3, Lb, 0, 1],
    [La, Lb, 0, 1]
])


patchTag = 1
P = 1
Q = 3

# Create a BSpline surface instance
surf = NURBS.Surface()

# Set surface degrees
surf.degree_u = P
surf.degree_v = Q

# Setting control points for surface
surf.set_ctrlpts(controlPts.tolist(), 2, 4)


# Set knot vectors
uKnot = generateKnotVector(surf.degree_u, surf.ctrlpts_size_u)
vKnot = generateKnotVector(surf.degree_v, surf.ctrlpts_size_v)

surf.knotvector_u = uKnot
surf.knotvector_v = vKnot

# Refine surface
operations.refine_knotvector(surf, [1, 2])
# Set surface degrees
surf.degree_u = 2
surf.degree_v = 5
# Set knot vectors refined
uKnot = generateKnotVector(surf.degree_u, surf.ctrlpts_size_u)
vKnot = generateKnotVector(surf.degree_v, surf.ctrlpts_size_v)
surf.knotvector_u = uKnot
surf.knotvector_v = vKnot

noPtsX = surf.ctrlpts_size_u
noPtsY = surf.ctrlpts_size_v

operations.insert_knot(surf,[None,0.95],[0,1])

# Visualize surface
surfVisualize(surf, hold=True)

# nDMaterial ElasticIsotropic $nDtag_elastic $elasticidad_probeta
# $poisson_probeta
E1 = 2.1e11  # Young's modulus N/m^2
E2 = E1
nu = 0.0  # Poisson's ratio
rho = 2.0e2  # kg/m^3


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

# matTags = [3, 4, 3, 4, 3]
# thickness = [10. * mm, 10. * mm, 10. * mm, 10. * mm, 10. * mm]
# θ = [0 * deg2rad, 45 * deg2rad, 90 * deg2rad, -45 * deg2rad, 0 * deg2rad]

matTags = [3]
thickness = [50. * mm]
θ = [0 * deg2rad]

gFact = [0.0, 0.0, 0.0]


Nlayers = len(θ)

controlPts = surf.ctrlpts2d[:] # Given in v,u
controlPts = np.array(compatibility.flip_ctrlpts2d(controlPts)) #Flipping tp u,v


print("controlPts.tolist(): ", controlPts.tolist())


ops.IGA("Patch", patchTag, surf.degree_u, surf.degree_v, surf.ctrlpts_size_u, surf.ctrlpts_size_v,
        "-type", "KLShell",
        # "-nonLinearGeometry", 0,
        "-planeStressMatTags", *matTags,
        "-gFact", *gFact,
        "-theta", *θ,
        "-thickness", *thickness,
        "-uKnot", *surf.knotvector_u, "-vKnot", *surf.knotvector_v, "-controlPts", *controlPts.flatten())

# exit()
nPoints = surf.ctrlpts_size_u * surf.ctrlpts_size_v

fixedNodes = np.concatenate((np.arange(1, surf.ctrlpts_size_u+1, 1), np.arange(surf.ctrlpts_size_u+1, surf.ctrlpts_size_v, 1)))
for n in ops.getNodeTags():
    if n in fixedNodes:
        ops.fix(n, 1, 1, 1)
    else:
        ops.fix(n, 0, 1, 0)

equalDOFnodes_master = np.arange(2*surf.ctrlpts_size_u+1, nPoints, surf.ctrlpts_size_u)
for masterNode in equalDOFnodes_master:
    retainedNodes = np.arange(masterNode+1, masterNode+surf.ctrlpts_size_u, 1)
    for retainedNode in retainedNodes:
        print("masterNode,retainedNode: ", masterNode,retainedNode)
        ops.equalDOF(int(masterNode), int(retainedNode), 1, 2, 3)




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
Pz = (2.7489e6)/4.0
followerLoadsPos = [0.0, 0.0, Pz/2] 
followerLoadsNeg = [0.0, 0.0, -Pz/2]

ops.eleLoad("-ele", 1, "-type", "-IGAFollowerLoad", 0.0, 1.0, *followerLoadsPos)
ops.eleLoad("-ele", 1, "-type", "-IGAFollowerLoad", 1.0, 1.0, *followerLoadsPos)

ops.eleLoad("-ele", 1, "-type", "-IGAFollowerLoad", 0.0, 0.95, *followerLoadsNeg)
ops.eleLoad("-ele", 1, "-type", "-IGAFollowerLoad", 1.0, 0.95, *followerLoadsNeg)
print("Finished loading nodes")

# ops.load(nPoints-2,0,0,Pz/2)
# ops.load(nPoints,0,0,Pz/2)



print("Starting analysis")

# create SOE
ops.system("FullGeneral")

# create DOF number
ops.numberer("Plain")

# create constraint handler
ops.constraints("Plain")
# ops.constraints("Transformation")

# create integrator
nSteps = 50
ops.integrator("LoadControl", 1.0 / nSteps)
# ops.integrator("LoadControl", 1.0)

# ops.algorithm("Linear")
# ops.algorithm("Newton")
# ops.algorithm("SecantNewton")
ops.algorithm("LineSearch")
# ops.algorithm("ModifiedNewton")
# ops.algorithm("KrylovNewton")

# Create test
ops.test("NormDispIncr", 1.0e-5, 300, 1)

# create analysis object
ops.analysis("Static")


# perform the analysis
import matplotlib.pyplot as plt
data = np.zeros((nSteps + 1, 2))
for j in range(nSteps):
    result=ops.analyze(1)
    if result!=0:
        break
        exit(-1)
    else:
        # Adding deformation to controlPts
        controlPts = surf.ctrlpts2d[:]
        controlPts = np.array(compatibility.flip_ctrlpts2d(controlPts)) #Flipping tp u,v

        fDef = 1e0
        i = 1
        for dim in controlPts:
            for point in dim:
                point[:3] += fDef * np.array(ops.nodeDisp(i))
                i += 1

        # Setting control points for surface
        controlPts = compatibility.flip_ctrlpts2d(controlPts)
        controlPts = (np.array(controlPts).reshape(
            nPoints, 4))
        surf.set_ctrlpts(controlPts.tolist(), surf.ctrlpts_size_u, surf.ctrlpts_size_v)

        # Visualize surface
        surfVisualize(surf, hold=True)


print("Done")




print("Finished analysis")


