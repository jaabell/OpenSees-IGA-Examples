#  IGA CANTILEVER PLATE UNDER SELF WEIGHT. ONE ELEMENT MESH, LINEAR CONVERGENCE OBTAINED


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
    [0      , 0  , 0 , 1] ,
    [La*1/4 , 0  , 0 , 1] ,
    [La*2/4 , 0  , 0 , 1] ,
    [La*3/4 , 0  , 0 , 1] ,
    [La     , 0  , 0 , 1] ,
    [0      , Lb , 0 , 1] ,
    [La*1/4 , Lb , 0 , 1] ,
    [La*2/4 , Lb , 0 , 1] ,
    [La*3/4 , Lb , 0 , 1] ,
    [La     , Lb , 0 , 1]
])


patchTag = 1
P = 1
Q = 4

# Create a BSpline surface instance
surf = NURBS.Surface()

# Set surface degrees
surf.degree_u = P
surf.degree_v = Q

# Setting control points for surface
surf.set_ctrlpts(controlPts.tolist(), 2, 5)

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
E1 = 2.1e11  # Young's modulus N/m^2
E2 = E1
nu = 0  # Poisson's ratio
rho = 1.0e4  # kg/m^3



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


gFact = [0.0, 0.0, 0.0]




Nlayers = len(θ)

controlPts = surf.ctrlpts2d[:]

print("controlPts =",controlPts)

controlPts = np.array(compatibility.flip_ctrlpts2d(controlPts))

print("controlPts =",controlPts)

print("controlPts.tolist(): ", controlPts.tolist())

print(uKnot)
print(vKnot)


nodeStartTag=1

ops.IGA("Patch", patchTag, nodeStartTag, P, Q, noPtsX, noPtsY,
        "-type", "KLShell",
        "-nonLinearGeometry", 0,
        "-planeStressMatTags", *matTags,
        "-gFact", *gFact,
        "-theta", *θ,
        "-thickness", *thickness,
        "-uKnot", *uKnot, "-vKnot", *vKnot, "-controlPts", *controlPts.flatten())



for n in [1,2,3,4,5,6,7,8,9,10]:
    if n in [1,2,3,4]:
        ops.fix(n,1,1,1)
    else:
        ops.fix(n,0,1,0)


print("\n\n\nPRINTING DOMAIN-----------------------")
ops.printModel()
print("\n\n\nDONE PRINTING DOMAIN-----------------------")

# ------------------------------
# Start of analysis generation
# ------------------------------

# create TimeSeries
# ops.timeSeries("Linear", 1)
ops.timeSeries("Constant", 1)

# create a plain load pattern
ops.pattern("Plain", 1, 1)

print("Loading nodes")
weight = [0.0, 0.0, -9.8066]
ops.eleLoad("-ele", 1, "-type", "-SelfWeight", *weight)
print("Finished loading nodes")


print("Calculating Rayleigh damping")
Nmodes = 2 # or whatever
w2 = ops.eigen(Nmodes)

# Pick your modes and damping ratios
wi = w2[0]**0.5; zetai = 0.05 # 2% in mode 1
wj = w2[1]**0.5; zetaj = 0.018 # 1.8% in mode 6

A = np.array([[1/wi, wi],[1/wj, wj]])
b = np.array([zetai,zetaj])

x = np.linalg.solve(A,2*b)

#             M    KT  KI  Kn
# ops.rayleigh(x[0],0.0,0.0,x[1])
ops.rayleigh(x[0],x[1],0.0,0.0)



print("Starting analysis")

# create SOE
ops.system("FullGeneral")

# create DOF number
ops.numberer("Plain")

# create constraint handler
ops.constraints("Plain")

# create integrator
ops.integrator("Newmark", 0.5, 0.25)

ops.algorithm("Linear")

# Create test
ops.test("NormDispIncr", 1.0e-8, 300,1)

# create analysis object
# ops.analysis("Static")
ops.analysis("Transient")

nSteps=1000
tMax=100

ops.recorder('Node', '-file', 'Node10.out', '-closeOnWrite','-time', '-node', *[10], '-dof', *[3], *['disp'])

# perform the analysis
import matplotlib.pyplot as plt
data=np.zeros((nSteps+1,2))
for j in range(nSteps):
    if ops.analyze(1,tMax/nSteps)!=0:
        print("analysis failed")

print("Finished analysis")




controlPts = surf.ctrlpts2d[:]
controlPts = compatibility.flip_ctrlpts2d(controlPts)

fDef = 1e0
i = 1
for dim in controlPts:
    for point in dim:
        point[:3] += fDef * np.array(ops.nodeDisp(i))
        i += 1

controlPts = compatibility.flip_ctrlpts2d(controlPts)

controlPts = (np.array(controlPts).reshape(
    surf.ctrlpts_size_u * surf.ctrlpts_size_v, 4)).tolist()

# Setting control points for surface
surf.set_ctrlpts(controlPts, 2, 5)


# Visualize surface
surfVisualize(surf, hold=True)

print("ops.nodeDisp(7,2): ", 1000*ops.nodeDisp(9,3),"mm")
print("ops.nodeDisp(8,2): ", 1000*ops.nodeDisp(10,3),"mm")

I=(Lb*(sum(thickness)**3))/12.0
W=rho*weight[2]*(sum(thickness)*Lb)
elasticSolution = abs(W*La**4/(8*E1*I))

print("elasticSolution: ", 1000*elasticSolution, "mm")

print("Done")
