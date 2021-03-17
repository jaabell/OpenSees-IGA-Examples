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
controlPtsL = surfL.ctrlpts2d[:]  # This gives v,u format
controlPtsL = np.array(compatibility.flip_ctrlpts2d(controlPtsL))  # This flips it to u,v
interfacePoints = controlPtsL[-1]
nodesOnLeftPatch = controlPtsL[-2]

controlPtsR = surfR.ctrlpts2d[:]  # This gives v,u format
controlPtsR = np.array(compatibility.flip_ctrlpts2d(controlPtsR))  # This flips it to u,v
nodesOnRightPatch = controlPtsR[1]

print("nodesOnLeftPatch: ", nodesOnLeftPatch)
print("interfacePoints: ", interfacePoints)
print("nodesOnRightPatch: ", nodesOnRightPatch)

# These are in u,v
controlPtsBendingStrip = np.concatenate([nodesOnLeftPatch, interfacePoints, nodesOnRightPatch])
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
ops.nDMaterial("ElasticOrthotropicPlaneStress", tagPlaneStress3, E1 * 30000000, E2 * 30000000, 0, 0, rho)
# ops.nDMaterial("ElasticOrthotropicPlaneStress", tagPlaneStress3, 0, 0, nu12, nu21, rho)


deg2rad = pi / 180

# matTags = [3, 4, 3, 4, 3]
# thickness = [10. * mm, 10. * mm, 10. * mm, 10. * mm, 10. * mm]
# θ = [0 * deg2rad, 45 * deg2rad, 90 * deg2rad, -45 * deg2rad, 0 * deg2rad]

matTags = [1]
thickness = [50. * mm]
θ = [0 * deg2rad]

gFact = [0.0, 0.0, 0.0 * 9.807]


Nlayers = len(θ)

names = ["leftPatch", "rightPatch", "bendingStrip"]
shellType = 'KLShell'

for i in range(len(container)):
    surf = container[i]
    name = names[i]
    if name == "bendingStrip":
        shellType = "KLShell_BendingStrip"

    # Flipping control point to u,v format
    controlPts = surf.ctrlpts2d[:]
    controlPts = np.array(compatibility.flip_ctrlpts2d(controlPts))

    # Creating a Patch in OpenSees
    ops.IGA("Patch", patchTag, nodeStartTag, surf.degree_u, surf.degree_v, surf.ctrlpts_size_u, surf.ctrlpts_size_v,
            "-type", shellType,
            "-nonLinearGeometry", 0,
            "-planeStressMatTags", *matTags,
            "-gFact", *gFact,
            "-theta", *θ,
            "-thickness", *thickness,
            "-uKnot", *surf.knotvector_u, "-vKnot", *surf.knotvector_v, "-controlPts", *controlPts.flatten())

    # Update patchTag, nodeStartTag and materialTag
    lastElTag = ops.getEleTags()[-1]
    lastNodeTag = ops.getNodeTags()[-1]
    matTags[0] += 1
    patchTag = lastElTag + 1
    nodeStartTag = lastNodeTag + 1


# Creating constraints
for n in ops.getNodeTags():
    n = int(n)
    if n in [1, 2, 3, 4]:
        ops.fix(n, 1, 1, 1)
    # else:
    #     ops.fix(n, 1, 1, 0)

# Creating equal dof's in patch interface
for n in [7, 8]:
    ops.equalDOF(n, n + 2, 1, 2, 3)

nodesOnPatches = [5, 6, 7, 8, 11, 12]
nodesOnBendingStrip = [17, 20, 18, 21, 19, 22]

for i in range(len(nodesOnPatches)):
    const = int(nodesOnPatches[i])
    ret = int(nodesOnBendingStrip[i])
    ops.equalDOF(const, ret, 1, 2, 3)


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
    ops.load(n, 0, 0, -Pz / 2.0)
    # ops.load(n, 0, Pz / 2.0, 0)
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


I = (Lb * (sum(thickness)**3)) / 12.0
elasticSolution = (Pz * ((2 * La)**3)) / (3 * E1 * I)

result = np.array(ops.nodeDisp(16, 3))
# result = np.array(ops.nodeDisp(8, 3))

print("result: ", result)
print("elasticSolution: ", elasticSolution)

nodesOnLeftPatch = [1, 2, 3, 4, 5, 6, 7, 8]
nodesOnRightPatch = [9, 10, 11, 12, 13, 14, 15, 16]
nodesOnBendingStrip = [17, 18, 19, 20, 21, 22]

controlPts_leftPatch = np.array(compatibility.flip_ctrlpts2d(surfL.ctrlpts2d[:])).reshape([surfL.ctrlpts_size_v * surfL.ctrlpts_size_u, 4]).tolist()
controlPts_rightPatch = np.array(compatibility.flip_ctrlpts2d(surfR.ctrlpts2d[:])).reshape([surfR.ctrlpts_size_v * surfR.ctrlpts_size_u, 4]).tolist()
controlPts_bendingStrip = np.array(compatibility.flip_ctrlpts2d(bendingStrip.ctrlpts2d[:])).reshape([bendingStrip.ctrlpts_size_v * bendingStrip.ctrlpts_size_u, 4]).tolist()

fdef = 1e2
for n in ops.getNodeTags():
    # Adding deformation to controlPts
    if n in nodesOnLeftPatch:
        indexN = nodesOnLeftPatch.index(n)
        point = controlPts_leftPatch[indexN]

    elif n in nodesOnRightPatch:
        indexN = nodesOnRightPatch.index(n)
        point = controlPts_rightPatch[indexN]

    elif n in nodesOnBendingStrip:
        indexN = nodesOnBendingStrip.index(n)
        point = controlPts_bendingStrip[indexN]

    point[0] = ops.nodeCoord(n)[0] + fdef * ops.nodeDisp(n)[0]
    point[1] = ops.nodeCoord(n)[1] + fdef * ops.nodeDisp(n)[1]
    point[2] = ops.nodeCoord(n)[2] + fdef * ops.nodeDisp(n)[2]


# # Setting control points for surfaces
# controlPts_leftPatch = compatibility.flip_ctrlpts2d(controlPts_leftPatch, surfL.ctrlpts_size_u, surf.ctrlpts_size_v)
# controlPts_rightPatch = compatibility.flip_ctrlpts2d(controlPts_rightPatch, surfR.ctrlpts_size_u, surfR.ctrlpts_size_v)
# controlPts_bendingStrip = compatibility.flip_ctrlpts2d(controlPts_bendingStrip, bendingStrip.ctrlpts_size_u, bendingStrip.ctrlpts_size_v)




surfL.degree_u=3
surfL.degree_v=1
surfL.set_ctrlpts(controlPts_leftPatch, 4, 2)
surfL.knotvector_u=knotvector.generate(surfL.degree_u, surfL.ctrlpts_size_u)
surfL.knotvector_v=knotvector.generate(surfL.degree_v, surfL.ctrlpts_size_v)

surfR.degree_u=3
surfR.degree_v=1
surfR.set_ctrlpts(controlPts_rightPatch, 4, 2)
surfR.knotvector_u=knotvector.generate(surfR.degree_u, surfR.ctrlpts_size_u)
surfR.knotvector_v=knotvector.generate(surfR.degree_v, surfR.ctrlpts_size_v)

bendingStrip.degree_u=1
bendingStrip.degree_v=2
bendingStrip.set_ctrlpts(controlPts_bendingStrip, 2, 3)
bendingStrip.knotvector_u=knotvector.generate(bendingStrip.degree_u, bendingStrip.ctrlpts_size_u)
bendingStrip.knotvector_v=knotvector.generate(bendingStrip.degree_v, bendingStrip.ctrlpts_size_v)

container.sample_size = 5
for surf in container:
    surf.evaluate()
# Visualization configuration
container.vis = VisMPL.VisSurface(ctrlpts=True, legend=False)
# container.vis.ctrlpts_offset=0.1

# Render the surface
evalcolor = ["green", "green", 'blue']
container.render(evalcolor=evalcolor)

