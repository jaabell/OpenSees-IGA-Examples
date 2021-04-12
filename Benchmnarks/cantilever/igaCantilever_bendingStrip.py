#  IGA CANTILEVER PLATE with bending strip


import numpy as np
import opensees as ops
from math import *

# Geomgl utilities for visualization and surface manipulation
from geomdl import NURBS, compatibility, operations, multi, knotvector
from surfVisualize import *

La = 10.0  	#
Lb = 2.0  	#
mm = 1.0 / 1000.  # m

ops.wipe()
ops.model('basic', '-ndm', 3, '-ndf', 3)

# These are given in v,u
controlPts = np.array([
    [0, 0, 0, 1],
    [0, Lb, 0, 1],
    [La * 1 / 3, 0, 0, 1],
    [La * 1 / 3, Lb, 0, 1],
    [La * 2 / 3, 0, 0, 1],
    [La * 2 / 3, Lb, 0, 1],
    [La, 0, 0, 1],
    [La, Lb, 0, 1]
])



# Create a NURBS surface instance
surfL = NURBS.Surface()

# Set surface degrees
surfL.degree_u = 3
surfL.degree_v = 1

# Setting control points for surface
surfL.set_ctrlpts(controlPts.tolist(), 4, 2)

# Set knot vectors
surfL.knotvector_u = knotvector.generate(surfL.degree_u, surfL.ctrlpts_size_u)
surfL.knotvector_v = knotvector.generate(surfL.degree_v, surfL.ctrlpts_size_v)


surfR = operations.translate(surfL, [La, 0, 0])


# Creating bending strip

from edgeHandler import *

interfacePoints, nodesOnLeftPatch=edgeGetter(surfL,"10","11")

nodesOnRightPatch = edgeGetter(surfR,"00","01")[1]

bendingStrip=makeBendingStrip(nodesOnLeftPatch,interfacePoints,nodesOnRightPatch)



# Creating container for multipatches

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


tagPlaneStress1 = 1
ops.nDMaterial("ElasticOrthotropicPlaneStress", tagPlaneStress1, E1, E2, nu12, nu21, rho)
# ops.nDMaterial("PlaneStress", tagPlaneStress1, tagNDmat1)

tagPlaneStress2 = 2
ops.nDMaterial("ElasticOrthotropicPlaneStress", tagPlaneStress2, E1, E2, nu12, nu21, rho)


tagPlaneStress3 = 3
ops.nDMaterial("ElasticOrthotropicPlaneStress", tagPlaneStress3, E1 * 1000, E2 * 1000, 0, 0, rho)


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
patchTag = 1
nodeStartTag = 1

nodesMap=[]

for i in range(len(container)):
    surf = container[i]
    name = names[i]

    # Identifying bending strip
    if name == "bendingStrip":
        shellType = "KLShell_BendingStrip"

    # Flipping control point to u,v format
    controlPts = surf.ctrlpts2d[:]
    controlPts = np.array(compatibility.flip_ctrlpts2d(controlPts))

    for dim in controlPts: # Unweighting control pts
        for point in dim:
            point[0:3]/=point[3]

    # Creating a Patch in OpenSees
    ops.IGA("Patch", patchTag, nodeStartTag, surf.degree_u, surf.degree_v, surf.ctrlpts_size_u, surf.ctrlpts_size_v,
            "-type", shellType,
            "-nonLinearGeometry", 0,
            "-planeStressMatTags", *matTags,
            "-gFact", *gFact,
            "-theta", *θ,
            "-thickness", *thickness,
            "-uKnot", *surf.knotvector_u, "-vKnot", *surf.knotvector_v, "-controlPts", *controlPts.flatten())


    # Get the nodes on current patch 
    nodesMap.append(np.arange(nodeStartTag,ops.getNodeTags()[-1]+1).tolist())

    # Update patchTag, nodeStartTag and materialTag
    lastElTag = ops.getEleTags()[-1]
    lastNodeTag = ops.getNodeTags()[-1]
    matTags[0] += 1
    patchTag = lastElTag + 1
    nodeStartTag = lastNodeTag + 1



# Creating constraints
for n in ops.getNodeTags():
    n = int(n)
    if n in [1, 2, 5, 6]:
        ops.fix(n, 1, 1, 1)


# Creating equal dof's in patch interface
for n in [4, 8]:
    ops.equalDOF(n, n + 5, 1, 2, 3)

nodesOnPatches = [3, 4, 10, 7, 8, 14]
nodesOnBendingStrip = [17, 18, 19, 20, 21, 22]

for i in range(len(nodesOnPatches)):
    const = int(nodesOnPatches[i])
    ret = int(nodesOnBendingStrip[i])
    ops.equalDOF(const, ret, 1, 2, 3)


print("ops.getEleTags(): ", ops.getEleTags())
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
# Cargar nodos 12,16
Pz = 1e2
for n in [12, 16]:
    # # Vertical load
    ops.load(n, 0, 0, -Pz / 2.0)
    # # Axial load
    # ops.load(n, -Pz/2.0, 0, 0)
print("Finished loading nodes")


print("Starting analysis")

# create SOE
ops.system("FullGeneral")

# create DOF number
ops.numberer("Plain")

# create constraint handler
ops.constraints("Plain")

# create integrator
ops.integrator("LoadControl", 1.0)

ops.algorithm("Linear")

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

print("result: ", result)
print("elasticSolution: ", elasticSolution)


# Adding deformation to control points
fdef = 1e2
for i in range(len(container)):
    surf=container[i]
    nodes=nodesMap[i]
    controlPts = np.array(compatibility.flip_ctrlpts2d(surf.ctrlpts2d[:])).reshape([surf.ctrlpts_size_v * surf.ctrlpts_size_u, 4]).tolist()
    for n in nodes:
        # Get node position in list
        indexN = nodes.index(n)
        point = controlPts[indexN]

        # Add deformation scaled by fdef
        weight=point[3]
        point[0] = (ops.nodeCoord(n)[0] + fdef * ops.nodeDisp(n)[0])*weight
        point[1] = (ops.nodeCoord(n)[1] + fdef * ops.nodeDisp(n)[1])*weight
        point[2] = (ops.nodeCoord(n)[2] + fdef * ops.nodeDisp(n)[2])*weight


    nPoints=surf.ctrlpts_size_u*surf.ctrlpts_size_v
    shape=np.array(compatibility.flip_ctrlpts2d(surf.ctrlpts2d[:])).shape
    controlPts=np.array(controlPts).reshape(shape)
    controlPts=np.array(compatibility.flip_ctrlpts2d(controlPts))

    surf.set_ctrlpts(controlPts.reshape(nPoints,4).tolist(), surf.ctrlpts_size_u, surf.ctrlpts_size_v)

container.sample_size = 10
for surf in container:
    surf.evaluate()
    
# Visualization configuration
container.vis = VisMPL.VisSurface(ctrlpts=True, legend=False, animate=True)
# container.vis = VisMPL.VisSurfWireframe(ctrlpts=True, legend=False)
# container.vis.ctrlpts_offset=0.1

# Render the surface
evalcolor = ["green", "green", 'black']
cpcolor=["green","green","black"]
container.render(evalcolor=evalcolor, cpcolor=cpcolor)

