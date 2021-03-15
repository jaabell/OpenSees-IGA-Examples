#  IGA CANTILEVER PLATE. GEOMETRY OBTAINED FROM NGUYEN'S igaCicleBendingStrip2D

# Validado, Lunes 1 de Marzo


import numpy as np
import opensees as ops
from math import *

# Geomgl utilities for visualization and surface manipulation
from geomdl import NURBS, compatibility, operations, knotvector
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


mm = 1.0 / 1000.  # m


ops.wipe()
ops.model('basic', '-ndm', 3, '-ndf', 3)


# These are given in v,u
controlPts = np.array([
    [4.953, 0., 0., 1.],
    [4.953, 0., 0.16172, 1.],
    [4.953, 0., 0.48516, 1.],
    [4.953, 0., 0.97031, 1.],
    [4.953, 0., 1.61719, 1.],
    [4.953, 0., 2.26406, 1.],
    [4.953, 0., 2.91094, 1.],
    [4.953, 0., 3.55781, 1.],
    [4.953, 0., 4.20469, 1.],
    [4.953, 0., 4.68984, 1.],
    [4.953, 0., 5.01328, 1.],
    [4.953, 0., 5.175, 1.],
    [4.953, 0.22298, 0., 0.98169],
    [4.953, 0.22298, 0.16172, 0.98169],
    [4.953, 0.22298, 0.48516, 0.98169],
    [4.953, 0.22298, 0.97031, 0.98169],
    [4.953, 0.22298, 1.61719, 0.98169],
    [4.953, 0.22298, 2.26406, 0.98169],
    [4.953, 0.22298, 2.91094, 0.98169],
    [4.953, 0.22298, 3.55781, 0.98169],
    [4.953, 0.22298, 4.20469, 0.98169],
    [4.953, 0.22298, 4.68984, 0.98169],
    [4.953, 0.22298, 5.01328, 0.98169],
    [4.953, 0.22298, 5.175, 0.98169],
    [4.92579, 0.68133, 0., 0.94813],
    [4.92579, 0.68133, 0.16172, 0.94813],
    [4.92579, 0.68133, 0.48516, 0.94813],
    [4.92579, 0.68133, 0.97031, 0.94813],
    [4.92579, 0.68133, 1.61719, 0.94813],
    [4.92579, 0.68133, 2.26406, 0.94813],
    [4.92579, 0.68133, 2.91094, 0.94813],
    [4.92579, 0.68133, 3.55781, 0.94813],
    [4.92579, 0.68133, 4.20469, 0.94813],
    [4.92579, 0.68133, 4.68984, 0.94813],
    [4.92579, 0.68133, 5.01328, 0.94813],
    [4.92579, 0.68133, 5.175, 0.94813],
    [4.79656, 1.38332, 0., 0.90695],
    [4.79656, 1.38332, 0.16172, 0.90695],
    [4.79656, 1.38332, 0.48516, 0.90695],
    [4.79656, 1.38332, 0.97031, 0.90695],
    [4.79656, 1.38332, 1.61719, 0.90695],
    [4.79656, 1.38332, 2.26406, 0.90695],
    [4.79656, 1.38332, 2.91094, 0.90695],
    [4.79656, 1.38332, 3.55781, 0.90695],
    [4.79656, 1.38332, 4.20469, 0.90695],
    [4.79656, 1.38332, 4.68984, 0.90695],
    [4.79656, 1.38332, 5.01328, 0.90695],
    [4.79656, 1.38332, 5.175, 0.90695],
    [4.4343, 2.3002, 0., 0.87033],
    [4.4343, 2.3002, 0.16172, 0.87033],
    [4.4343, 2.3002, 0.48516, 0.87033],
    [4.4343, 2.3002, 0.97031, 0.87033],
    [4.4343, 2.3002, 1.61719, 0.87033],
    [4.4343, 2.3002, 2.26406, 0.87033],
    [4.4343, 2.3002, 2.91094, 0.87033],
    [4.4343, 2.3002, 3.55781, 0.87033],
    [4.4343, 2.3002, 4.20469, 0.87033],
    [4.4343, 2.3002, 4.68984, 0.87033],
    [4.4343, 2.3002, 5.01328, 0.87033],
    [4.4343, 2.3002, 5.175, 0.87033],
    [3.87817, 3.15152, 0., 0.85203],
    [3.87817, 3.15152, 0.16172, 0.85203],
    [3.87817, 3.15152, 0.48516, 0.85203],
    [3.87817, 3.15152, 0.97031, 0.85203],
    [3.87817, 3.15152, 1.61719, 0.85203],
    [3.87817, 3.15152, 2.26406, 0.85203],
    [3.87817, 3.15152, 2.91094, 0.85203],
    [3.87817, 3.15152, 3.55781, 0.85203],
    [3.87817, 3.15152, 4.20469, 0.85203],
    [3.87817, 3.15152, 4.68984, 0.85203],
    [3.87817, 3.15152, 5.01328, 0.85203],
    [3.87817, 3.15152, 5.175, 0.85203],
    [3.15152, 3.87817, 0., 0.85203],
    [3.15152, 3.87817, 0.16172, 0.85203],
    [3.15152, 3.87817, 0.48516, 0.85203],
    [3.15152, 3.87817, 0.97031, 0.85203],
    [3.15152, 3.87817, 1.61719, 0.85203],
    [3.15152, 3.87817, 2.26406, 0.85203],
    [3.15152, 3.87817, 2.91094, 0.85203],
    [3.15152, 3.87817, 3.55781, 0.85203],
    [3.15152, 3.87817, 4.20469, 0.85203],
    [3.15152, 3.87817, 4.68984, 0.85203],
    [3.15152, 3.87817, 5.01328, 0.85203],
    [3.15152, 3.87817, 5.175, 0.85203],
    [2.3002, 4.4343, 0., 0.87033],
    [2.3002, 4.4343, 0.16172, 0.87033],
    [2.3002, 4.4343, 0.48516, 0.87033],
    [2.3002, 4.4343, 0.97031, 0.87033],
    [2.3002, 4.4343, 1.61719, 0.87033],
    [2.3002, 4.4343, 2.26406, 0.87033],
    [2.3002, 4.4343, 2.91094, 0.87033],
    [2.3002, 4.4343, 3.55781, 0.87033],
    [2.3002, 4.4343, 4.20469, 0.87033],
    [2.3002, 4.4343, 4.68984, 0.87033],
    [2.3002, 4.4343, 5.01328, 0.87033],
    [2.3002, 4.4343, 5.175, 0.87033],
    [1.38332, 4.79656, 0., 0.90695],
    [1.38332, 4.79656, 0.16172, 0.90695],
    [1.38332, 4.79656, 0.48516, 0.90695],
    [1.38332, 4.79656, 0.97031, 0.90695],
    [1.38332, 4.79656, 1.61719, 0.90695],
    [1.38332, 4.79656, 2.26406, 0.90695],
    [1.38332, 4.79656, 2.91094, 0.90695],
    [1.38332, 4.79656, 3.55781, 0.90695],
    [1.38332, 4.79656, 4.20469, 0.90695],
    [1.38332, 4.79656, 4.68984, 0.90695],
    [1.38332, 4.79656, 5.01328, 0.90695],
    [1.38332, 4.79656, 5.175, 0.90695],
    [0.68133, 4.92579, 0., 0.94813],
    [0.68133, 4.92579, 0.16172, 0.94813],
    [0.68133, 4.92579, 0.48516, 0.94813],
    [0.68133, 4.92579, 0.97031, 0.94813],
    [0.68133, 4.92579, 1.61719, 0.94813],
    [0.68133, 4.92579, 2.26406, 0.94813],
    [0.68133, 4.92579, 2.91094, 0.94813],
    [0.68133, 4.92579, 3.55781, 0.94813],
    [0.68133, 4.92579, 4.20469, 0.94813],
    [0.68133, 4.92579, 4.68984, 0.94813],
    [0.68133, 4.92579, 5.01328, 0.94813],
    [0.68133, 4.92579, 5.175, 0.94813],
    [0.22298, 4.953, 0., 0.98169],
    [0.22298, 4.953, 0.16172, 0.98169],
    [0.22298, 4.953, 0.48516, 0.98169],
    [0.22298, 4.953, 0.97031, 0.98169],
    [0.22298, 4.953, 1.61719, 0.98169],
    [0.22298, 4.953, 2.26406, 0.98169],
    [0.22298, 4.953, 2.91094, 0.98169],
    [0.22298, 4.953, 3.55781, 0.98169],
    [0.22298, 4.953, 4.20469, 0.98169],
    [0.22298, 4.953, 4.68984, 0.98169],
    [0.22298, 4.953, 5.01328, 0.98169],
    [0.22298, 4.953, 5.175, 0.98169],
    [0., 4.953, 0., 1.],
    [0., 4.953, 0.16172, 1.],
    [0., 4.953, 0.48516, 1.],
    [0., 4.953, 0.97031, 1.],
    [0., 4.953, 1.61719, 1.],
    [0., 4.953, 2.26406, 1.],
    [0., 4.953, 2.91094, 1.],
    [0., 4.953, 3.55781, 1.],
    [0., 4.953, 4.20469, 1.],
    [0., 4.953, 4.68984, 1.],
    [0., 4.953, 5.01328, 1.],
    [0., 4.953, 5.175, 1.]
])

for point in controlPts:  # weighting points (x*w,y*w,z*w,w)
    for i in range(3):
        point[i] *= point[3]

patchTag = 1
P = 5
Q = 5

# Create a BSpline surface instance
surf = NURBS.Surface()

# Set surface degrees
surf.degree_u = P
surf.degree_v = Q

# Setting control points for surface
surf.set_ctrlpts(controlPts.tolist(), 12, 12)

# Set knot vectors
surf.knotvector_u = knotvector.generate(surf.degree_u, surf.ctrlpts_size_u)
surf.knotvector_v = knotvector.generate(surf.degree_v, surf.ctrlpts_size_v)


# Visualize surface
surfVisualize(surf, hold=True)


# nDMaterial ElasticIsotropic $nDtag_elastic $elasticidad_probeta
# $poisson_probeta
E1 = 1.05e7  # Young's modulus N/m^2
E2 = E1
nu = 0.3125  # Poisson's ratio
rho = 0.5e2  # kg/m^3
t = 0.094


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
thickness = [t]
θ = [0 * deg2rad]

gFact = [0.0, 0.0, 0.0]


Nlayers = len(θ)

controlPts = surf.ctrlpts2d[:]  # Given in v,u
controlPts = np.array(compatibility.flip_ctrlpts2d(
    controlPts))  # Flipping to u,v

for point in controlPts:
    for i in range(surf.ctrlpts_size_v):
        for k in range(3):
            point[i][k] /= point[i][3]


ops.IGA("Patch", patchTag, surf.degree_u, surf.degree_v, surf.ctrlpts_size_u, surf.ctrlpts_size_v,
        "-type", "KLShell",
        # "-nonLinearGeometry", 0,
        "-planeStressMatTags", *matTags,
        "-gFact", *gFact,
        "-theta", *θ,
        "-thickness", *thickness,
        "-uKnot", *surf.knotvector_u, "-vKnot", *surf.knotvector_v, "-controlPts", *controlPts.flatten())


nPoints = surf.ctrlpts_size_u * surf.ctrlpts_size_v
print("surf.ctrlpts_size_u: ", surf.ctrlpts_size_u)
print("surf.ctrlpts_size_v: ", surf.ctrlpts_size_v)

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



for n in ops.getNodeTags():
    if n in nodesOnAD:
        ops.fix(n, 1, 0, 0)
    elif n in nodesOnCB:
        ops.fix(n, 0, 1, 0)
    elif n in nodesOnDC:
        ops.fix(n, 0, 0, 1)


for i in range(len(nodesOnAD)):
    # ops.equalDOF(int(nodesOnAD[i]), int(nodesNextToAD[i]), 1, 2, 3)
    ops.equalDOF(int(nodesOnAD[i]), int(nodesNextToAD[i]), 2)
    # ops.equalDOF(int(nodesOnCB[i]), int(nodesNextToCB[i]), 1, 2, 3)
    ops.equalDOF(int(nodesOnCB[i]), int(nodesNextToCB[i]), 2)
    # ops.equalDOF(int(nodesOnDC[i]), int(nodesNextToDC[i]), 1, 2, 3)
    ops.equalDOF(int(nodesOnDC[i]), int(nodesNextToDC[i]), 3)


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

Pz = 40000.0
forcedNode = surf.ctrlpts_size_u
ops.load(forcedNode, 0, Pz / 4.0, 0)


print("Finished loading nodes")


print("Starting analysis")

# Create test
ops.test("NormDispIncr", 1.0e-5, 200, 1) # Apparently faster
# ops.test("NormUnbalance", 1.0e-4, 80, 1)
# ops.test("EnergyIncr", 1.0e-4, 80, 1)

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
# ops.algorithm("NewtonLineSearch", 'type', 'Bisection')
ops.algorithm("NewtonLineSearch")
# ops.algorithm("ModifiedNewton")
# ops.algorithm("KrylovNewton")
# ops.algorithm("BFGS")
# ops.algorithm("Broyden")

# create integrator
# delta = 0.1
# defMax = 2.7
# nSteps = abs(int(defMax / delta))

nSteps = 40
ops.integrator("LoadControl", 1.0 / nSteps)
# ops.integrator("DisplacementControl", forcedNode, 2, delta)

# create analysis object
ops.analysis("Static")


# perform the analysis
import matplotlib.pyplot as plt
data = np.zeros((nSteps + 1, 4))

for j in range(nSteps):
    print("=================================")
    print(f"Load step {j}")
    print("=================================")
    result = ops.analyze(1)
    if result != 0:
        break
        exit(-1)
    else:
        # Adding deformation to controlPts
        controlPts = surf.ctrlpts2d[:]
        controlPts = compatibility.flip_ctrlpts2d(
            controlPts)  # Flipping to u,v

        fDef = 1
        i = 1
        for dim in controlPts:
            for point in dim:
                # Times the weight
                point[:3] += fDef * np.array(ops.nodeDisp(i)) * point[3]
                i += 1

        # Setting control points for surface
        controlPts = compatibility.flip_ctrlpts2d(controlPts)
        controlPts = (np.array(controlPts).reshape(
            nPoints, 4))
        surf.set_ctrlpts(controlPts.tolist(),
                         surf.ctrlpts_size_u, surf.ctrlpts_size_v)

        # Visualize surface
        surfVisualize(surf, hold=True)

        controlPts = surf.ctrlpts2d[:]
        controlPts = compatibility.flip_ctrlpts2d(
            controlPts)  # Flipping to u,v
        i = 1
        for dim in controlPts:
            for point in dim:
                point[:3] -= fDef * np.array(ops.nodeDisp(i)) * point[3]
                i += 1

        # Setting control points for surface
        controlPts = compatibility.flip_ctrlpts2d(controlPts)
        controlPts = (np.array(controlPts).reshape(
            nPoints, 4))
        surf.set_ctrlpts(controlPts.tolist(),
                         surf.ctrlpts_size_u, surf.ctrlpts_size_v)

        data[j + 1, 0] = abs(ops.nodeDisp(forcedNode, 2))
        data[j + 1, 1] = ops.getLoadFactor(1)
        data[j + 1, 2] = abs(ops.nodeDisp(1, 1))  # Point B
        data[j + 1, 3] = abs(ops.nodeDisp(133, 1))  # Point C
        # elasticSolution = (data[j + 1, 2] * (La**3)) / (3 * E1 * I)
        # print("ops.getLoadFactor(1)*Pz: ", ops.getLoadFactor(1) * Pz)
        # print("elasticSolution: ", elasticSolution)
        print("data[j+1,0]: ", data[j + 1, 0])
        print("data[j+1,1]: ", data[j + 1, 1])
        print("4*data[j+1,1]: ", 4 * data[j + 1, 1])

        # B=ops.printB('-ret')
        # print("B: ", B)

        print("\nNext load step\n")

# % refenrece solution from Sze et al, 2004
# % Popular benchmark problems for geometric nonlinear analysis of shells

P = Pz/1000 * np.array([0, 0.025 , 0.05  , 0.075 , 0.1   , 0.15  , 0.20  , 0.25  , 0.30  , 0.35  , 0.40  , 0.45  , 0.5   , 0.525 , 0.55  , 0.6   , 0.65  , 0.7   , 0.75  , 0.8   , 0.85  , 0.9   , 0.95  , 1])
wA =     np.array([0, 0.819      , 1.26  , 1.527 , 1.707 , 1.936 , 2.079 , 2.180 , 2.257 , 2.321 , 2.376 , 2.425 , 2.473 , 2.543 , 2.577 , 2.618 , 2.648 , 2.672 , 2.692 , 2.710 , 2.726 , 2.741 , 2.755 , 2.768])
uB =     np.array([0, 0.864      , 1.471 , 1.901 , 2.217 , 2.641 , 2.904 , 3.087 , 3.227 , 3.342 , 3.443 , 3.539 , 3.653 , 4.061 , 4.171 , 4.274 , 4.338 , 4.385 , 4.423 , 4.455 , 4.483 , 4.508 , 4.53  , 4.551])
uC =     np.array([0, 0.872      , 1.493 , 1.946 , 2.293 , 2.792 , 3.106 , 3.310 , 3.452 , 3.556 , 3.632 , 3.688 , 3.718 , 3.580 , 3.518 , 3.452 , 3.410 , 3.378 , 3.353 , 3.332 , 3.313 , 3.297 , 3.283 , 3.269])


plt.plot(data[:, 0], Pz /1000* data[:, 1], 'or')
plt.plot(data[:, 2], Pz /1000* data[:, 1], 'ob')
plt.plot(data[:, 3], Pz /1000* data[:, 1], 'og')


plt.plot(wA, P, '-r')
plt.plot(uB, P, '-b')
plt.plot(uC, P, '-g')

plt.xlabel('Upward deflection at point A')
plt.ylabel('Pulling force at point A')
plt.show()

# Visualize surface
# surfVisualize(surf, hold=True)

# print("Done")

# elasticSolution = (Pz * (La**3)) / (3 * E1 * I)

# print("elasticSolution: ", elasticSolution)
# print("data[nSteps,0]: ", data[nSteps, 0])


# print("Finished analysis")
