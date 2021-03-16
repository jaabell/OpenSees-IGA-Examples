#  IGA CANTILEVER PLATE with bending strip


import numpy as np
import opensees as ops
from math import *

# Geomgl utilities for visualization and surface manipulation
from geomdl import NURBS, compatibility, operations, multi, knotvector
from surfVisualize import *

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
nodeStartTag = 1
P = 1
Q = 3

# Create a NURBS surface instance
surfL = NURBS.Surface()

# Set surface degrees
surfL.degree_u = P
surfL.degree_v = Q

# Setting control points for surface
surfL.set_ctrlpts(controlPts.tolist(), 2, 4)

# Set knot vectors
uKnot = knotvector.generate(surfL.degree_u, surfL.ctrlpts_size_u)
vKnot = knotvector.generate(surfL.degree_v, surfL.ctrlpts_size_v)

surfL.knotvector_u = uKnot
surfL.knotvector_v = vKnot


noPtsX = surfL.ctrlpts_size_u
noPtsY = surfL.ctrlpts_size_v

surfR = operations.translate(surfL, [La, 0, 0])


# Creating bending strip
controlPtsL = surfL.ctrlpts2d[:] # This gives v,u format
controlPtsL = np.array(compatibility.flip_ctrlpts2d(controlPtsL)) # This flips it to u,v
interfacePoints= controlPtsL[-1]
nodesOnLeftPatch = controlPtsL[-2]

controlPtsR = surfR.ctrlpts2d[:] # This gives v,u format
controlPtsR = np.array(compatibility.flip_ctrlpts2d(controlPtsR)) # This flips it to u,v
nodesOnRightPatch = controlPtsR[1]

print("nodesOnLeftPatch: ", nodesOnLeftPatch)
print("interfacePoints: ", interfacePoints)
print("nodesOnRightPatch: ", nodesOnRightPatch)

# These are in u,v
controlPtsBendingStrip=np.concatenate([nodesOnLeftPatch,interfacePoints,nodesOnRightPatch])
print("controlPtsBendingStrip: ", controlPtsBendingStrip)

bendingStrip = NURBS.Surface()

# Set surface degrees
bendingStrip.degree_u = 2
bendingStrip.degree_v = 1

# Setting control points for surface
bendingStrip.set_ctrlpts(controlPtsBendingStrip.tolist(), 3, 2)

# Set knot vectors
bendingStrip.knotvector_u = knotvector.generate(bendingStrip.degree_u, bendingStrip.ctrlpts_size_u)
bendingStrip.knotvector_v = knotvector.generate(bendingStrip.degree_v, bendingStrip.ctrlpts_size_v)

surfList = [surfL, surfR, bendingStrip]

container = multi.SurfaceContainer(surfList)

# Visualize surface

container.sample_size = 5
for surf in container:
    surf.evaluate()
# Visualization configuration
container.vis = VisMPL.VisSurface(ctrlpts=True, legend=False)
# container.vis.ctrlpts_offset=0.1

# Render the surface
evalcolor = ["green", "green", 'blue']
container.render(evalcolor=evalcolor)

# exit()
# surfVisualize(surfList, hold=True)


# nDMaterial ElasticIsotropic $nDtag_elastic $elasticidad_probeta
# $poisson_probeta
E1 = 2.1e11  # Young's modulus N/m^2
E2 = E1
nu12 = 0.0  # Poisson's ratio
nu21 = 0.0  # Poisson's ratio
rho = 8.0e3  # *9.807 # kg/m^3
t = 0.05


# nDMaterial PlateFiber $nDtag_platefiber $nDtag_elastic
tagPlaneStress1 = 1
# ops.nDMaterial("ElasticPlaneStress", tagPlaneStress1, E1, nu, rho)
ops.nDMaterial("ElasticOrthotropicPlaneStress", tagPlaneStress1, E1, E2, nu12, nu21, rho)
# ops.nDMaterial("PlaneStress", tagPlaneStress1, tagNDmat1)

tagPlaneStress2 = 2
# ops.nDMaterial("PlaneStress", tagPlaneStress2, tagNDmat2)
ops.nDMaterial("ElasticOrthotropicPlaneStress", tagPlaneStress2, E1, E2, nu12, nu21, rho)


tagPlaneStress3 = 3
# ops.nDMaterial("PlaneStress", tagPlaneStress2, tagNDmat2)
ops.nDMaterial("ElasticOrthotropicPlaneStress", tagPlaneStress3, E1, E2, nu12, nu21, rho)


deg2rad = pi / 180

# matTags = [3, 4, 3, 4, 3]
# thickness = [10. * mm, 10. * mm, 10. * mm, 10. * mm, 10. * mm]
# θ = [0 * deg2rad, 45 * deg2rad, 90 * deg2rad, -45 * deg2rad, 0 * deg2rad]

matTags = [1]
thickness = [50. * mm]
θ = [0 * deg2rad]

gFact = [0.0, 0.0, 0.0 * 9.807]


Nlayers = len(θ)

names=["leftPatch","rightPatch","bendingStrip"]

for i in range(len(container)):
    surf=container[i]
    name=names[i]
    shellType='KLShell'
    if name=="bendingStrip":
        shellType="KLShell_BendingStrip"

    # Flipping control point to u,v format
    controlPts = surf.ctrlpts2d[:]
    controlPts = np.array(compatibility.flip_ctrlpts2d(controlPts))

    # Creating a Patch in OpenSees
    ops.IGA("Patch", patchTag, nodeStartTag, P, Q, noPtsX, noPtsY,
            "-type", shellType,
            "-nonLinearGeometry", 0,
            "-planeStressMatTags", *matTags,
            "-gFact", *gFact,
            "-theta", *θ,
            "-thickness", *thickness,
            "-uKnot", *uKnot, "-vKnot", *vKnot, "-controlPts", *controlPts.flatten())

    # Update patchTag, nodeStartTag and materialTag
    lastElTag = ops.getEleTags()[-1]
    lastNodeTag = ops.getNodeTags()[-1]
    matTags[0] += 1
    patchTag = lastElTag + 1
    nodeStartTag = lastNodeTag + 1


# Creating constraints
for n in np.arange(1, 17):
    n = int(n)
    if n in [1, 2, 3, 4]:
        ops.fix(n, 1, 1, 1)
    else:
        ops.fix(n, 1, 1, 0)

# Creating equal dof's in patch interface
for n in [9, 10]:
    ops.equalDOF(n, n + 2, 1, 2, 3)


print("ops.getEleTags(): ", ops.getEleTags())
print("\n\n\nPRINTING DOMAIN-----------------------")
ops.printModel()
print("\n\n\nDONE PRINTING DOMAIN-----------------------")
# exit()

# ------------------------------
# Start of analysis generation
# ------------------------------

# create TimeSeries
ops.timeSeries("Linear", 1)

# create a plain load pattern
ops.pattern("Plain", 1, 1)

print("Loading nodes")
# Cargar nodos 7,8
Pz = 1e2
# for n in [7, 8]:
for n in [15, 16]:
    # ops.load(n,0,0,-Pz/2.0)
    ops.load(n, 0, Pz / 2.0, 0)
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
# ops.integrator("LoadControl", 1.0)

ops.algorithm("Linear")
# ops.algorithm("Newton")
# ops.algorithm("ModifiedNewton")
# ops.algorithm("KrylovNewton")

# Create test
ops.test("NormDispIncr", 1.0e-9, 50, 1)

# create analysis object
ops.analysis("Static")


# perform the analysis
ops.analyze(1)


print("Finished analysis!")


# Haven't coded deformation for multipatch yet
exit()

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
surf.set_ctrlpts(controlPts, 2, 4)


# Set knot vectors
uKnot = generateKnotVector(surf.degree_u, surf.ctrlpts_size_u)
vKnot = generateKnotVector(surf.degree_v, surf.ctrlpts_size_v)

surf.knotvector_u = uKnot
surf.knotvector_v = vKnot

noPtsX = surf.ctrlpts_size_u
noPtsY = surf.ctrlpts_size_v

# Visualize surface
surfVisualize(surf, hold=True)

print("ops.nodeDisp(7,2): ", 1000 * np.array(ops.nodeDisp(7)), "mm")
print("ops.nodeDisp(8,2): ", 1000 * np.array(ops.nodeDisp(8)), "mm")

I = (Lb * (sum(thickness)**3)) / 12.0
elasticSolution = (Pz * (La**3)) / (3 * E1 * I)

result = 1000 * np.array(ops.nodeDisp(7, 3))

print("result: ", result)
print("elasticSolution: ", 1000 * elasticSolution)

print("Done")
