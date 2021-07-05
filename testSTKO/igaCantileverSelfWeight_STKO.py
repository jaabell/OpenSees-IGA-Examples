#  IGA CANTILEVER PLATE UNDER SELF WEIGHT. ONE ELEMENT MESH, LINEAR CONVERGENCE OBTAINED

import numpy as np
import opensees as ops
from math import *


La = 10.0   # v length
Lb = 1.0    # u length
mm = 1.0 / 1000.  # meter to milimeter

ops.wipe()
ops.model('basic', '-ndm', 3, '-ndf', 3)  # 3D model, 3 DOF per node


# These are in u, v order
controlPts = np.array([
    [0, 0, 0, 1],
    [0, Lb, 0, 1],
    [La*1/4, 0, 0, 1],
    [La*1/4, Lb, 0, 1],
    [La*2/4, 0, 0, 1],
    [La*2/4, Lb, 0, 1],
    [La*3/4, 0, 0, 1],
    [La*3/4, Lb, 0, 1],
    [La, 0, 0, 1],
    [La, Lb, 0, 1]
])


# Material parameters
E1 = 2.1e11  # Young's modulus N/m^2
E2 = E1
nu = 0  # Poisson's ratio
rho = 1.0e4  # kg/m^3


# Creating the necessary materials
tagNDmat1 = 10
ops.nDMaterial("ElasticIsotropic", tagNDmat1, E1, nu, rho)

tagNDmat2 = 20
ops.nDMaterial("ElasticIsotropic", tagNDmat2, E2, nu, rho)


# # Creatring the plane stress used for the element formulation
# tagPlaneStress1 = 3
# ops.nDMaterial("PlaneStress", tagPlaneStress1, tagNDmat1)

# tagPlaneStress2 = 4
# ops.nDMaterial("PlaneStress", tagPlaneStress2, tagNDmat2)

# Material parameters
MPa = 1e6
E1 = 2.1e11  # Young's modulus N/m^2
E2 = E1
nu12 = 0  # Poisson's ratio
nu21 = 0
G12 = E1/(2*(1+nu12))
rho = 1.0e4  # kg/m^3

Xt = 350.0*MPa
Xc = 350.0*MPa
Yt = 300.0*MPa
Yc = 300.0*MPa
S = 90.0*MPa
c1 = 4.0e-6
c2 = 30.0
c3 = 2.0e-6
c4 = 0.8
c5 = 80.0
c6 = 0.0
c7 = 0.0
c8 = 0.0
c9 = 0.0
b = 1.0

tagPlaneStress1 = 3
vonPaepeParams = [E1, E2, nu12, nu21, G12, rho, 
            Xt, Xc, Yt, Yc, S, c1, c2, c3, c4, c5, c6, c7, c8, c9, b]
ops.nDMaterial("VonPapaDamage", tagPlaneStress1, *vonPaepeParams)

tagPlaneStress2 = 4
vonPaepeParams = [E1, E2, nu12, nu21, G12, rho, 
            Xt, Xc, Yt, Yc, S, c1, c2, c3, c4, c5, c6, c7, c8, c9, b]
ops.nDMaterial("VonPapaDamage", tagPlaneStress2, *vonPaepeParams)


deg2rad = pi / 180  # for conversion from radian to degrees

matTags = [tagPlaneStress1, tagPlaneStress2, tagPlaneStress1, tagPlaneStress2, tagPlaneStress1]  # material tag for each layer
thickness = [10. * mm, 10. * mm, 10. * mm, 10. *
             mm, 10. * mm]  # Thickness of each layer
θ = [0 * deg2rad, 45 * deg2rad, 90 * deg2rad, -45 * deg2rad,
     0 * deg2rad]  # Angle of orientation of each layer

# body force acceleration factors (optional, using selfWeight here)
gFact = [0.0, 0.0, 0.0]

# knot vectors along each direction
uKnot = [0, 0, 1, 1]
vKnot = [0, 0, 0, 0, 0, 1, 1, 1, 1, 1]

# NURBS basis functions order (P in "u" direction, Q in "v" direction)
P = 1
Q = 4

# Number of control points in each direction
noPtsX = 2
noPtsY = 5

# Number of layers/plies for each IGA patch
Nlayers = len(θ)

# The patch is an element, so the element tags start with tag "patchTag + 1". Intended for multiPatches
patchTag = 1

# The tag of the first node/control point in the patch (in case of multiPatches, have to update this with the last added node)
nodeStartTag = 1

# secType = 'Elastic'
# secTag = 1
# secArgs = [E1, 10, 100]
# ops.section('PlateFiber', 1, 10, 0.1)

# ops.IGA call to create the patch (super element)
ops.IGA("Patch", patchTag, nodeStartTag, P, Q, noPtsX, noPtsY,
        "-type", "KLShell", # Element type to use for the patch (used when creating a bending strip)
        # Flag that tells whether to use linear or non-linear geometry, linear in this case
        "-nonLinearGeometry", 1,
        "-planeStressMatTags", *matTags,
        "-gFact", *gFact,
        "-theta", *θ,
        "-thickness", *thickness,
        "-uKnot", *uKnot, "-vKnot", *vKnot, "-controlPts", *controlPts.flatten())



# Boundary conditions, fixing first two rows for "clamping" condition, fixing "y" displacement for the rest
for n in ops.getNodeTags():
    if n in [1, 2, 3, 4]:
        ops.fix(n, 1, 1, 1)
    else:
        ops.fix(n, 0, 1, 0)


print("\n\n\nPRINTING DOMAIN-----------------------")
ops.printModel()
print("\n\n\nDONE PRINTING DOMAIN-----------------------")

# STKO Recorder

ops.recorder("Element", "-xml", "fiberstuff.xml","-time","-ele",2,"material","1","fiber","0","stress")
ops.recorder("mpco","iga_cantilever",
    "-N","displacement",
    "-E","section.fiber.stress",
    "-E","section.fiber.strain",
    "-E","section.fiber.damagestate",
    )

print(f"DONE! ")

# exit(0)

# ------------------------------
# Start of analysis generation
# ------------------------------

# create TimeSeries
ops.timeSeries("Linear", 1)

# create a plain load pattern
ops.pattern("Plain", 1, 1)

# Loading the patch, selfweight in this case
weight = [0.0, 0.0, -9.8066]
# -ele 1 means that the patch with tag 1 will load it's associated elements
ops.eleLoad("-ele", 1, "-type", "-SelfWeight", *weight)


print("Starting analysis")

# create SOE
ops.system("FullGeneral")

# create DOF number
ops.numberer("Plain")

# create constraint handler
ops.constraints("Plain")

# create integrator
ops.integrator("LoadControl", 1.0)

# create algorithm
ops.test("NormDispIncr", 1.0e-8, 50, 2)
ops.algorithm("Newton")

# create analysis object
ops.analysis("Static")


# perform the analysis
ops.analyze(1)

print("Finished analysis!\n")

print("Measured vertical displacements at the tip : ")
print("ops.nodeDisp(7,2): ", 1000*ops.nodeDisp(9, 3), "mm")
print("ops.nodeDisp(8,2): ", 1000*ops.nodeDisp(10, 3), "mm\n")


I = (Lb*(sum(thickness)**3))/12.0
W = rho*weight[2]*(sum(thickness)*Lb)
elasticSolution = W*La**4/(8*E1*I)

ops.record()

print("elasticSolution: ", 1000*elasticSolution, "mm\n")

print("Done analysis!")

ops.remove("recorders")

