#  IGA CANTILEVER ROTATING


import numpy as np
import opensees as ops
from math import *

# Geomgl utilities for visualization and surface manipulation
from geomdl import NURBS, compatibility, operations, knotvector
from surfVisualize import *


La = 6    #
Lb = 0.5
mm = 1.0 / 1000.  # m

d = 0.2  # radius

ops.wipe()
ops.model('basic', '-ndm', 3, '-ndf', 3)

# These are given in v,u
controlPts = np.array([
    [0.00000 * La/10 + d, 0.0 * Lb, 0.0, 1.0],
    [0.00000 * La/10 + d, 1.0 * Lb, 0.0, 1.0],
    [0.31250 * La/10 + d, 0.0 * Lb, 0.0, 1.0],
    [0.31250 * La/10 + d, 1.0 * Lb, 0.0, 1.0],
    [0.93750 * La/10 + d, 0.0 * Lb, 0.0, 1.0],
    [0.93750 * La/10 + d, 1.0 * Lb, 0.0, 1.0],
    [1.87500 * La/10 + d, 0.0 * Lb, 0.0, 1.0],
    [1.87500 * La/10 + d, 1.0 * Lb, 0.0, 1.0],
    [3.12500 * La/10 + d, 0.0 * Lb, 0.0, 1.0],
    [3.12500 * La/10 + d, 1.0 * Lb, 0.0, 1.0],
    [4.37500 * La/10 + d, 0.0 * Lb, 0.0, 1.0],
    [4.37500 * La/10 + d, 1.0 * Lb, 0.0, 1.0],
    [5.62500 * La/10 + d, 0.0 * Lb, 0.0, 1.0],
    [5.62500 * La/10 + d, 1.0 * Lb, 0.0, 1.0],
    [6.87500 * La/10 + d, 0.0 * Lb, 0.0, 1.0],
    [6.87500 * La/10 + d, 1.0 * Lb, 0.0, 1.0],
    [8.12500 * La/10 + d, 0.0 * Lb, 0.0, 1.0],
    [8.12500 * La/10 + d, 1.0 * Lb, 0.0, 1.0],
    [9.06250 * La/10 + d, 0.0 * Lb, 0.0, 1.0],
    [9.06250 * La/10 + d, 1.0 * Lb, 0.0, 1.0],
    [9.68750 * La/10 + d, 0.0 * Lb, 0.0, 1.0],
    [9.68750 * La/10 + d, 1.0 * Lb, 0.0, 1.0],
    [10.00000*La/10 + d, 0.0 * Lb, 0.0, 1.0],
    [10.00000*La/10 + d, 1.0 * Lb, 0.0, 1.0]
])

# controlPts[0][0]*=0
# controlPts[1][0]*=0
# controlPts[2][0]*=0
# controlPts[3][0]*=0

# controlPts[0][2]*=0
# controlPts[1][2]*=0
# controlPts[2][2]*=0
# controlPts[3][2]*=0

patchTag = 1
nodeStartTag = 1
P = 5
Q = 1

# Create a BSpline surface instance
surf = NURBS.Surface()

# Set surface degrees
surf.degree_u = P
surf.degree_v = Q

# Setting control points for surface
surf.set_ctrlpts(controlPts.tolist(), 12, 2)

# Set knot vectors
uKnot = knotvector.generate(surf.degree_u, surf.ctrlpts_size_u)
vKnot = knotvector.generate(surf.degree_v, surf.ctrlpts_size_v)

surf.knotvector_u = uKnot
surf.knotvector_v = vKnot


noPtsX = surf.ctrlpts_size_u
noPtsY = surf.ctrlpts_size_v

# Visualize surface
surfVisualize(surf, hold=True)


# nDMaterial ElasticIsotropic $nDtag_elastic $elasticidad_probeta
# $poisson_probeta
GPa = 1e9
E1 = 5*GPa  # Young's modulus N/m^2
E2 = E1
nu = 0.0  # Poisson's ratio
rho = 8e3  # *9.807 # kg/m^3
t = 0.5


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


matTags = [3,4]
thickness = [t/2, t/2]
θ = [0 * deg2rad] * 2
g = 9.807

I = (Lb*(sum(thickness)**3))/12.0
W = rho*g*(sum(thickness)*Lb)
elasticSolution = abs(W*La**4/(8*E1*I))


gFact = [0.0, 0.0, 0*g]


Nlayers = len(θ)

controlPts = surf.ctrlpts2d[:]
controlPts = np.array(compatibility.flip_ctrlpts2d(controlPts))


ops.IGA("Patch", patchTag, nodeStartTag, P, Q, noPtsX, noPtsY,
        "-type", "KLShell",
        # "-nonLinearGeometry", 0,
        "-planeStressMatTags", *matTags,
        "-gFact", *gFact,
        "-theta", *θ,
        "-thickness", *thickness,
        "-uKnot", *uKnot, "-vKnot", *vKnot, "-controlPts", *controlPts.flatten())

# exit()

# Creating constraints
for n in ops.getNodeTags():
    if n in [1, 2, 13, 14]:
        # ops.fix(n, 1, 1, 1)
        pass
    else:
        ops.fix(n, 0, 1, 0)

nSpins = 3
ω = 5  # rad/s
tMax = (nSpins*2*pi/ω)  # 5 seconds
deltaT = 0.05/4
t = np.arange(0, tMax+deltaT, deltaT)


# First row parametrization

dX_0 = d*np.cos(ω*t) 
dX_0 -= dX_0[0]
dZ_0 = d*np.sin(ω*t)
dZ_0 -= dZ_0[0]


# Second row parametrization
d_1 = d + 0.03125*La  # distance to middle
dX_1 = d_1*np.cos(ω*t)
dX_1 -= dX_1[0] 
dZ_1 = d_1*np.sin(ω*t)
dZ_1 -= dZ_1[0] 



# create TimeSeries for first row

# in X
ops.timeSeries("Path", 1, '-time', *(t.tolist()), '-values', *(dX_0.tolist()), '-prependZero')

# create a plain load pattern
ops.pattern("Plain", 1, 1)
# ops.pattern("UniformExcitation", 1, -'disp', 1, '-vel0', 0.0)

# creating sp constraints
ops.sp(1, 1, 1.0)
ops.sp(13, 1, 1.0)


# in Z
ops.timeSeries("Path", 2, '-time', *(t.tolist()), '-values', *(dZ_0.tolist()), '-prependZero')

# create a plain load pattern
ops.pattern("Plain", 2, 2)

# creating sp constraints
ops.sp(1, 3, 1.0)
ops.sp(13, 3, 1.0)


# create TimeSeries for second row

# in X
ops.timeSeries("Path", 3, '-time', *(t.tolist()), '-values', *(dX_1.tolist()), '-prependZero')

# create a plain load pattern
ops.pattern("Plain", 3, 3)

# creating sp constraints
ops.sp(2, 1, 1.0)
ops.sp(14, 1, 1.0)


# in Z
ops.timeSeries("Path", 4, '-time', *(t.tolist()), '-values', *(dZ_1.tolist()), '-prependZero')

# create a plain load pattern
ops.pattern("Plain", 4, 4)

# creating sp constraints
ops.sp(2, 3, 1.0)
ops.sp(14, 3, 1.0)


# print("\n\n\nPRINTING DOMAIN-----------------------")
ops.printModel()
# print("\n\n\nDONE PRINTING DOMAIN-----------------------")


# from matplotlib.pylab import *

# plot(t,dX_0)
# plot(t,dX_1)
# show()

# exit()


# ------------------------------
# Start of analysis generation
# ------------------------------

# create TimeSeries
ops.timeSeries("Constant", 5)


# create a plain load pattern
ops.pattern("Plain", 5, 5)

weight = [0.0, 0.0, -g]
ops.eleLoad("-ele", 1, "-type", "-SelfWeight", *weight)


# Analysis
ops.test("EnergyIncr", 1.0e-8, 80, 0)
# ops.test("NormUnbalance", 1.0e-10, 90, 0)
# ops.test("NormDispIncr", 1.0e-8, 50, 0)

ops.algorithm("Newton")
# ops.algorithm("Linear")
# ops.algorithm("NewtonLineSearch")
# ops.algorithm("NewtonLineSearch", 'type', 'Bisection')

ops.numberer("RCM")

# ops.constraints("Plain")
ops.constraints("Transformation")

ops.integrator("Newmark", 0.5, 0.25)
ops.system("UmfPack")
ops.analysis("Transient")

# Create recorder
ops.recorder('Node', '-file', 'Node12_Z.out', '-closeOnWrite', '-time', '-node', *[12], '-dof', *[3], *['disp'])
ops.recorder('Node', '-file', 'Node12_X.out', '-closeOnWrite', '-time', '-node', *[12], '-dof', *[1], *['disp'])
# ops.recorder('Node', '-file', 'Node.out', '-time', '-node', *[12], '-dof', *[3], *['disp'])

nSteps = int(tMax/deltaT)


for j in range(nSteps):
    print(j/nSteps*100, '%')
    if (ops.analyze(1, deltaT) != 0):
        exit()

    node=1
    print(np.array(ops.nodeCoord(1))+np.array(ops.nodeDisp(1)))

    if j % 10 == 0:

        controlPts = surf.ctrlpts2d[:]
        controlPts = compatibility.flip_ctrlpts2d(controlPts)

        i = 1
        fdef = 1e0
        for dim in controlPts:
            for point in dim:
                point[0] = ops.nodeCoord(i, 1) + fdef * ops.nodeDisp(i, 1)
                point[1] = ops.nodeCoord(i, 2) + fdef * ops.nodeDisp(i, 2)
                point[2] = ops.nodeCoord(i, 3) + fdef * ops.nodeDisp(i, 3)
                i += 1

        nPoints = surf.ctrlpts_size_u * surf.ctrlpts_size_v
        shape = np.array(compatibility.flip_ctrlpts2d(surf.ctrlpts2d[:])).shape
        controlPts = np.array(controlPts).reshape(shape)
        controlPts = np.array(compatibility.flip_ctrlpts2d(controlPts))

        # print(controlPts)

        surf.set_ctrlpts(controlPts.reshape(nPoints, 4).tolist(),
                         surf.ctrlpts_size_u, surf.ctrlpts_size_v)

        # Visualize surface
        surfVisualize(surf, hold=False)

print(ops.nodeDisp(8, 3))

# ops.record()
print("Done")
