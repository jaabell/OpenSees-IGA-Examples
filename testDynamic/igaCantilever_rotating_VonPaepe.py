#  IGA CANTILEVER ROTATING


import numpy as np
import opensees as ops
from math import *

# Geomgl utilities for visualization and surface manipulation
from geomdl import NURBS, compatibility, operations, knotvector
from surfVisualize import *


La = 2.26    #
Lb = 0.15


d = 0.2  # radius

ops.wipe()
ops.model('basic', '-ndm', 3, '-ndf', 3)

# These are given in v,u
controlPts = np.array([
    [0.00000 / 10 * La + d, 0.0 * Lb, 0.0, 1.0],
    [0.00000 / 10 * La + d, 1.0 * Lb, 0.0, 1.0],
    [0.31250 / 10 * La + d, 0.0 * Lb, 0.0, 1.0],
    [0.31250 / 10 * La + d, 1.0 * Lb, 0.0, 1.0],
    [0.93750 / 10 * La + d, 0.0 * Lb, 0.0, 1.0],
    [0.93750 / 10 * La + d, 1.0 * Lb, 0.0, 1.0],
    [1.87500 / 10 * La + d, 0.0 * Lb, 0.0, 1.0],
    [1.87500 / 10 * La + d, 1.0 * Lb, 0.0, 1.0],
    [3.12500 / 10 * La + d, 0.0 * Lb, 0.0, 1.0],
    [3.12500 / 10 * La + d, 1.0 * Lb, 0.0, 1.0],
    [4.37500 / 10 * La + d, 0.0 * Lb, 0.0, 1.0],
    [4.37500 / 10 * La + d, 1.0 * Lb, 0.0, 1.0],
    [5.62500 / 10 * La + d, 0.0 * Lb, 0.0, 1.0],
    [5.62500 / 10 * La + d, 1.0 * Lb, 0.0, 1.0],
    [6.87500 / 10 * La + d, 0.0 * Lb, 0.0, 1.0],
    [6.87500 / 10 * La + d, 1.0 * Lb, 0.0, 1.0],
    [8.12500 / 10 * La + d, 0.0 * Lb, 0.0, 1.0],
    [8.12500 / 10 * La + d, 1.0 * Lb, 0.0, 1.0],
    [9.06250 / 10 * La + d, 0.0 * Lb, 0.0, 1.0],
    [9.06250 / 10 * La + d, 1.0 * Lb, 0.0, 1.0],
    [9.68750 / 10 * La + d, 0.0 * Lb, 0.0, 1.0],
    [9.68750 / 10 * La + d, 1.0 * Lb, 0.0, 1.0],
    [10.00000 / 10 * La + d, 0.0 * Lb, 0.0, 1.0],
    [10.00000 / 10 * La + d, 1.0 * Lb, 0.0, 1.0]
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
P = 3
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

factor = 1
E1 = factor*56101.68868e6  # Young's modulus N/m^2
E2 = factor*28050.680868e6
nu12 = 0.055  # Poisson's ratio
nu21 = E2*nu12/E1  # Poisson's ratio
G12 = factor*4488.174327e6  # Shear modulus
rho = 3181.27  # *9.807 # kg/m^3

Xt = 70e6
Xc = 70e6
Yt = 70e6
Yc = 70e6
S = 70e6
c1 = 4.0e-6
c2 = 30
c3 = 2.0e-6
c4 = 0.8
c5 = 80
c6 = 0
c7 = 0
c8 = 0
c9 = 0
b = 1

tagPlaneStress = 1
vonPaepeParams = [E1, E2, nu12, nu21, G12, rho, Xt, Xc,
                  Yt, Yc, S, c1, c2, c3, c4, c5, c6, c7, c8, c9, b]
ops.nDMaterial("VonPapaDamage", tagPlaneStress, *vonPaepeParams)
print("Created Von-Paepe material")
deg2rad = pi / 180


thickness = 2*1e0*np.array([2.43e-4, 2.43e-4, 2.43e-4, 2.43e-4, 2.43e-4, 2.43e-4, 2.43e-4, 2.43e-4, 2.43e-4, 2.43e-4,
                            2.43e-4, 2.43e-4, 2.43e-4, 2.43e-4, 2.43e-4, 2.43e-4, 2.43e-4, 2.43e-4, 2.43e-4, 2.43e-4, 2.43e-4, 2.43e-4])
θ = 1*deg2rad*np.array([45.0, -45.0, 0.0, 90.0, 0.0, 90.0, 0.0, 90.0, 45.0, -45.0,
                        0.0, 90.0, 0.0, 90.0, 45.0, -45.0, 45.0, -45.0, 45.0, -45.0, 45.0, -45.0])
Nlayers = len(θ)
matTags = [1]*Nlayers


controlPts = surf.ctrlpts2d[:]
controlPts = np.array(compatibility.flip_ctrlpts2d(controlPts))


ops.IGA("Patch", patchTag, nodeStartTag, P, Q, noPtsX, noPtsY,
        "-type", "KLShell",
        # "-nonLinearGeometry", 0,
        "-planeStressMatTags", *matTags,
        "-theta", *θ,
        "-thickness", *thickness,
        "-uKnot", *uKnot, "-vKnot", *vKnot, "-controlPts", *controlPts.flatten())





# ops.setParameter('-val', 0, "-ele", 1, "resetMaxStress")

# ops.updateParameter()

print("Succesfully created model")
# exit()

# Creating constraints
for n in ops.getNodeTags():
    if n in [1, 2, 13, 14]:
        # ops.fix(n, 1, 1, 1)
        pass
    else:
        ops.fix(n, 0, 1, 0)

nSpins = 3
ω = 1.5  # rad/s
tMax = (nSpins*2*pi/ω)  # 5 seconds
deltaT = 1e-3
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

g = 9.807
weight = [0.0, 0.0, -g]
ops.eleLoad("-ele", 1, "-type", "-SelfWeight", *weight)


# Analysis
ops.test("EnergyIncr", 1.0e-7, 100, 0)
# ops.test("NormUnbalance", 1.0e-7, 100, 0)
# ops.test("NormDispIncr", 1.0e-7, 100, 0)

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
    nCycles = 10
    ops.setParameter('-val', int(nCycles), "-ele", 1, "advanceDamageState")

    node = 1
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
        # surfVisualize(surf, hold=False)

print(ops.nodeDisp(8, 3))

# ops.record()
print("Done")
