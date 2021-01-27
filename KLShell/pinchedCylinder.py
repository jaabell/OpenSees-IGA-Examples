from geomdl.shapes import surface
from geomdl import exchange, compatibility, operations, helpers

from surfVisualize import *

import numpy as np
import opensees as ops
from math import *


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

cm = 1 / 100  # m

# Generate cylindirical surface
surf = surface.cylinder(radius=9 * cm, height=30 * cm)

# Set evaluation delta
surf.delta = 0.01

# Evaluate the surface
surf.evaluate()


# Refine knot vectors and update geometry (refine both directions)
# u,v
refU = 0
refV = 0
operations.refine_knotvector(surf, [refU, refV])


# Visualize surface
surfVisualize(surf, hold=True)


# Obtain surface geometry
print("Obtaining surface geometry")
uKnot = surf.knotvector_u
vKnot = surf.knotvector_v
P = surf.degree_u
Q = surf.degree_v
noPtsX = surf.ctrlpts_size_u
noPtsY = surf.ctrlpts_size_v
weights = surf.weights
controlPts = surf.ctrlpts2d[:]
controlPts = np.array(compatibility.flip_ctrlpts2d(controlPts))
print(controlPts)


print("uKnot: ", uKnot)
print("vKnot: ", vKnot)
print("P: ", P)
print("Q: ", Q)

exit()

#
controlPts = [
[[0.09       , 0.        , 0.	 , 1.   ] , [0.09      , 0.        , 0.3     , 1.     ] ],
[[0.063639   , 0.063639  , 0.  , 0.7071 ] , [0.063639  , 0.063639  , 0.21213 , 0.7071 ] ],
[[0.         , 0.09      , 0.  , 1.     ] , [0.        , 0.09      , 0.3     , 1.     ] ],
[[-0.063639  , 0.063639  , 0.  , 0.7071 ] , [-0.063639 , 0.063639  , 0.21213 , 0.7071 ] ],
[[-0.09      , 0.        , 0.  , 1.     ] , [-0.09     , 0.        , 0.3     , 1.     ] ],
[[-0.063639  , -0.063639 , 0.  , 0.7071 ] , [-0.063639 , -0.063639 , 0.21213 , 0.7071 ] ],
[[0.       	 , -0.09     , 0.  , 1.     ] , [0.        , -0.09     , 0.3     , 1.     ] ],
[[0.063639 	 , -0.063639 , 0.  , 0.7071 ] , [0.063639  , -0.063639 , 0.21213 , 0.7071 ] ],
[[0.09       , 0.        , 0.  , 1.     ] , [0.09      , 0.        , 0.3     , 1.     ] ]
]

#

ops.wipe()
ops.model('basic', '-ndm', 3, '-ndf', 3)

# Defining material properties
E1 = 2.1e11  # Young's modulus N/m^2
E2 = E1
nu = 0.0  # Poisson's ratio
rho = 8.0e3  # *9.807 # kg/m^3

mm = 1.0 / 1000.  # m
deg2rad = pi / 180

# nDMaterial ElasticIsotropic $nDtag_elastic $elasticidad_probeta
tagNDmat1 = 1
ops.nDMaterial("ElasticIsotropic", tagNDmat1, E1, nu, rho)

tagNDmat2 = 2
ops.nDMaterial("ElasticIsotropic", tagNDmat2, E2, nu, rho)


# nDMaterial PlateFiber $nDtag_platefiber $nDtag_elastic
tagPlaneStress1 = 3
ops.nDMaterial("PlaneStress", tagPlaneStress1, tagNDmat1)

tagPlaneStress2 = 4
ops.nDMaterial("PlaneStress", tagPlaneStress2, tagNDmat2)


matTags = [3, 4, 3, 4, 3]
thickness = [10. * mm, 10. * mm, 10. * mm, 10. * mm, 10. * mm]
θ = [0 * deg2rad, 45 * deg2rad, 90 * deg2rad, -45 * deg2rad, 0 * deg2rad]


gFact = [0.0, 0.0, 0.0 * 9.807]


# Numbering of patch
patchTag = 1


ops.IGA("Patch", patchTag, P, Q, noPtsX, noPtsY,
        "-type", "KLShell",
        # "-nonLinearGeometry", 0,
        "-planeStressMatTags", *matTags,
        "-gFact", *gFact,
        "-theta", *θ,
        "-thickness", *thickness,
        "-uKnot", *uKnot, "-vKnot", *vKnot, "-controlPts", *controlPts.flatten())


print("\n\n\nPRINTING DOMAIN-----------------------")
ops.printModel()
print("\n\n\nDONE PRINTING DOMAIN-----------------------")

# Fijar nodos 1, 2 Apoyo base
for n in [1, 2]:
    ops.fix(n, 1, 1, 1)

# Fijar nodos 9, 10 Superior fijos en x, z
for n in [9, 10]:
    ops.fix(n, 0, 1, 1)

# Equal dof en base 1,2 con 17,18
ops.equalDOF(1, 17, 1, 2, 3)
ops.equalDOF(2, 18, 1, 2, 3)

# ------------------------------
# Start of analysis generation
# ------------------------------

# create TimeSeries
ops.timeSeries("Linear", 1)

# create a plain load pattern
ops.pattern("Plain", 1, 1)

print("Loading nodes")
# Cargar nodos 7,8
Px = 0.5e0
for n in [9, 10]:
    ops.load(n, Px, 0, 0)
print("Finished loading nodes")


print("Starting analysis")

# create SOE
ops.system("FullGeneral")

# create DOF number
ops.numberer("Plain")

# create constraint handler
ops.constraints("Plain")

# create integrator
nSteps = 1
ops.integrator("LoadControl", 1.0 / nSteps)

# Create algorithm
# ops.algorithm("Linear")
ops.algorithm("Newton")


# Create test
ops.test("NormDispIncr", 1.0e-9, 50, 1)

# create analysis object
ops.analysis("Static")


# perform the analysis
import matplotlib.pyplot as plt
data = np.zeros((nSteps + 1, 2))
for j in range(nSteps):
    ops.analyze(1)
    data[j + 1, 0] = abs(1000 * ops.nodeDisp(9, 1))
    data[j + 1, 1] = ops.getLoadFactor(1) * (2 * Px)
    print("data[j+1,0],data[j+1,1]: ", data[j + 1, 0], data[j + 1, 1])

plt.plot(data[:, 0], data[:, 1])
plt.plot(data[:, 0], data[:, 1], 'or')
plt.xlabel('Vertical Displacement')
plt.ylabel('Vertical Load')
plt.show()


# # Upadting geometry with deformation
# oldShape=controlPts.shape
# controlPts = np.reshape(controlPts,(noPtsY*noPtsX,4))
# nodes=ops.getNodeTags()
# for node in nodes:
# 	controlPts[node-1,:3] += np.array(ops.nodeDisp(node))
# controlPts = np.reshape(controlPts,oldShape)

# # surf.set_ctrlpts(controlPts, noPtsX, noPtsY)
# controlPts = (np.array(controlPts).reshape(
#     surf.ctrlpts_size_u * surf.ctrlpts_size_v, 4)).tolist()

# # Setting control points for surface
# surf.set_ctrlpts(controlPts, noPtsX, noPtsY)

controlPts = surf.ctrlpts2d[:]
controlPts = compatibility.flip_ctrlpts2d(controlPts)

fDef = 1e5
i = 1
for dim in controlPts:
    for point in dim:
        point[:3] += fDef * np.array(ops.nodeDisp(i))
        i += 1

controlPts = compatibility.flip_ctrlpts2d(controlPts)

controlPts = (np.array(controlPts).reshape(
    surf.ctrlpts_size_u * surf.ctrlpts_size_v, 4)).tolist()

# Setting control points for surface
surf.set_ctrlpts(controlPts, noPtsX, noPtsY)

# Visualize surface


print("controlPts: ", controlPts)

surfVisualize(surf, hold=True)


print("Finished analysis")