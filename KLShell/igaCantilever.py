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




La = 50.0  	#
Lb = 1.0  	#
mm = 1.0 / 1000.  # m

ops.wipe()
ops.model('basic', '-ndm', 3, '-ndf', 3)

# These are given in v,u
controlPts = np.array([ 
    [0      , 0  , 0 , 2] ,
    [La*1/3 , 0  , 0 , 2] ,
    [La*2/3 , 0  , 0 , 2] ,
    [La     , 0  , 0 , 2] ,
    [0      , Lb , 0 , 2] ,
    [La*1/3 , Lb , 0 , 2] ,
    [La*2/3 , Lb , 0 , 2] ,
    [La     , Lb , 0 , 2]
])

controlPts = np.array([ 
    [0      , 0  , 0 , 2] ,
    [0      , Lb , 0 , 2] ,
    [La*1/3 , 0  , 0 , 2] ,
    [La*1/3 , Lb , 0 , 2] ,
    [La*2/3 , 0  , 0 , 2] ,
    [La*2/3 , Lb , 0 , 2] ,
    [La     , 0  , 0 , 2] ,
    [La     , Lb , 0 , 2]
])


patchTag = 1
P = 3
Q = 1

# Create a BSpline surface instance
surf = NURBS.Surface()

# Set surface degrees
surf.degree_u = P
surf.degree_v = Q

# Setting control points for surface
surf.set_ctrlpts(controlPts.tolist(), 4, 2)

# Set knot vectors
uKnot = generateKnotVector(surf.degree_u, surf.ctrlpts_size_u)
vKnot = generateKnotVector(surf.degree_v, surf.ctrlpts_size_v)

surf.knotvector_u = uKnot
surf.knotvector_v = vKnot


noPtsX = surf.ctrlpts_size_u
noPtsY = surf.ctrlpts_size_v

# Visualize surface
surfVisualize(surf, hold=True)


# nDMaterial ElasticIsotropic $nDtag_elastic $elasticidad_probeta
# $poisson_probeta
E1 = 100  # Young's modulus N/m^2
E2 = E1
nu = 0.0  # Poisson's ratio
rho = 8.0e3  # *9.807 # kg/m^3
t = 1


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

matTags = [3, 4, 3, 4, 3]
thickness = [10. * mm, 10. * mm, 10. * mm, 10. * mm, 10. * mm]
θ = [0 * deg2rad, 45 * deg2rad, 90 * deg2rad, -45 * deg2rad, 0 * deg2rad]

matTags = [3]
thickness = [t]
θ = [0 * deg2rad]

gFact = [0.0, 0.0, 0.0 * 9.807]


Nlayers = len(θ)

controlPts = surf.ctrlpts2d[:]
controlPts = np.array(compatibility.flip_ctrlpts2d(controlPts))


ops.IGA("Patch", patchTag, P, Q, noPtsX, noPtsY,
        "-type", "KLShell",
        "-nonLinearGeometry", 0,
        "-planeStressMatTags", *matTags,
        "-gFact", *gFact,
        "-theta", *θ,
        "-thickness", *thickness,
        "-uKnot", *uKnot, "-vKnot", *vKnot, "-controlPts", *controlPts.flatten())

# exit()


for n in [1,2,3,4,5,6,7,8]:
    if n in [1,2,5,6]:
        ops.fix(n,1,1,1)
    else:
        ops.fix(n,1,1,0)

# # #Fijar nodos 1, 2, 3, 4
# for n in [1,2,3,4]:
#     ops.fix(n,1,1,1)

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
# Cargar nodos 7,8
Pz=1
for n in [4,8]:
    ops.load(n,0,0,-Pz)
print("Finished loading nodes")





print("Starting analysis")

# create SOE
ops.system("UmfPack")

# create DOF number
ops.numberer("Plain")

# create constraint handler
ops.constraints("Plain")

# create integrator
nSteps=1
ops.integrator("LoadControl", 1.0/nSteps)
# ops.integrator("LoadControl", 1.0)

ops.algorithm("Linear")
# ops.algorithm("Newton")
# ops.algorithm("ModifiedNewton")
# ops.algorithm("KrylovNewton")

# Create test
ops.test("NormDispIncr", 1.0e-9, 50,1)

# create analysis object
ops.analysis("Static")


# perform the analysis
import matplotlib.pyplot as plt
data=np.zeros((nSteps+1,2))
for j in range(nSteps):
    ops.analyze(1)
    data[j+1,0] = ops.nodeDisp(4,3)
    data[j+1,1] = ops.getLoadFactor(1)*(2*Pz)
    # print("data[j+1,0],data[j+1,1]: ", data[j+1,0],data[j+1,1])

plt.plot(data[:,0], data[:,1])
plt.plot(data[:,0], data[:,1],'or')
plt.xlabel('Vertical Displacement')
plt.ylabel('Vertical Load')
plt.show()

print("Finished analysis")

controlPts = surf.ctrlpts2d[:]
controlPts = compatibility.flip_ctrlpts2d(controlPts)

fDef = 1e2
i = 1
for dim in controlPts:
    for point in dim:
        point[:3] += fDef * np.array(ops.nodeDisp(i)) 
        i += 1

controlPts = compatibility.flip_ctrlpts2d(controlPts)

controlPts = (np.array(controlPts).reshape(
    surf.ctrlpts_size_u * surf.ctrlpts_size_v, 4)).tolist()

# Setting control points for surface
surf.set_ctrlpts(controlPts, 4, 2)


# Set knot vectors
uKnot = generateKnotVector(surf.degree_u, surf.ctrlpts_size_u)
vKnot = generateKnotVector(surf.degree_v, surf.ctrlpts_size_v)

surf.knotvector_u = uKnot
surf.knotvector_v = vKnot

noPtsX = surf.ctrlpts_size_u
noPtsY = surf.ctrlpts_size_v

# Visualize surface
surfVisualize(surf, hold=True)

print("ops.nodeDisp(7,2): ", np.array(ops.nodeDisp(4)),"mm")
print("ops.nodeDisp(8,2): ", np.array(ops.nodeDisp(8)),"mm")

I=(Lb*(sum(thickness)**3))/12.0
elasticSolution=(2*Pz*(La**3))/(3*E1*I)

result=np.array(ops.nodeDisp(4,3))

print("result: ", result)
print("elasticSolution: ", elasticSolution)

print("Done")

