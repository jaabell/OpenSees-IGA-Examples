
from scipy.sparse.linalg import eigsh
from geomdl.visualization import VisVTK, VisMPL
from edgeHandler import *
from scipy import interpolate
import numpy as np
import opensees as ops
from math import *

# For saving dictionaries using pickle

from six.moves import cPickle as pickle #for performance

def save_dict(di_, filename_):
    with open(filename_, 'wb') as f:
        pickle.dump(di_, f)



# Geomdl utilities for visualization and surface manipulation
from geomdl import NURBS, compatibility, operations, multi, knotvector


from fileGetter import surfFromFile

# For each surface, I need it's control points, order, knot vectors, startNodeTag, endNodeTag, type
surfaces_data = {}

surfList = []
topSurf_original = surfFromFile("coarse.1.dat")
bottomSurf_original = surfFromFile("coarse.2.dat")


# Interpolating length, in order to know where to split the surfaces

# Getting the bounding box, average between two surfaces
boundingBox = (np.array(topSurf_original.bbox) +
               np.array(bottomSurf_original.bbox))/2
zLimits = boundingBox[:, 2]  # Getting the limits for parametrization
zLength = sum(abs(zLimits))

parametricInt = interpolate.interp1d(zLimits, [0, 1])
uCuts = np.array([0.011, 0.251, 0.696, 1.141, 1.586]
                 )  # Places to define sections

# Parametric coordinates (along "u") that define the sections
parametricCuts = parametricInt(uCuts)

surfList = []

topSurf = topSurf_original
bottomSurf = bottomSurf_original

surfTag=0
nodeStartTag = 1

for i in range(len(uCuts)):
    if i == 0:
        topSurf = topSurf_original
        bottomSurf = bottomSurf_original

    # Updating the interpolator to know where to cut next
    # Getting the bounding box, average between two surfaces
    boundingBox = (np.array(topSurf.bbox)+np.array(bottomSurf.bbox))/2
    zLimits = boundingBox[:, 2]  # Getting the limits for parametrization
    zLength = sum(abs(zLimits))

    parametricInt = interpolate.interp1d(zLimits, [0, 1])

    uCut = float(parametricInt(uCuts[i]))  # Parametric coordinate to cut at
    # This returns two surfaces, [0,uCut] and [uCut,1]
    cutSurfaces_top = operations.split_surface_u(topSurf, param=uCut)
    # This returns two surfaces, [0,uCut] and [uCut,1]
    cutSurfaces_bottom = operations.split_surface_u(bottomSurf, param=uCut)

    # Adding cutted surfaces to surfList
    if i == range(len(uCuts))[-1]:
        # surf_top = cutSurfaces_top[1]
        # surf_bottom = cutSurfaces_bottom[1]

        surf_top = cutSurfaces_top[0]
        surf_bottom = cutSurfaces_bottom[0]
    else:
        surf_top = cutSurfaces_top[0]
        surf_bottom = cutSurfaces_bottom[0]

    surfList.append(surf_top)
    surfList.append(surf_bottom)


    # Now adding surfaces to the global surfaces dictionary

    surfaces=[surf_top,surf_bottom]

    for surf in surfaces:

        surf_data = {}

        # Flipping control point to u,v format
        surf_data["nodeStartTag"] = nodeStartTag
        surf_data["ctrlpts"] = np.array(
            compatibility.flip_ctrlpts2d(surf.ctrlpts2d[:]))
        surf_data["uKnot"] = surf.knotvector_u
        surf_data["vKnot"] = surf.knotvector_v
        surf_data["p"] = surf.degree_u
        surf_data["q"] = surf.degree_v
        surf_data["nPtsU"] = surf.ctrlpts_size_u
        surf_data["nPtsV"] = surf.ctrlpts_size_v
        surf_data["nodeEndTag"] = surf_data["nodeStartTag"] + \
            surf_data["nPtsU"]*surf_data["nPtsV"] - 1
        surf_data["type"] = "Patch"

        # Creating node mapping
        nodeEndTag= surf_data["nodeEndTag"]
        nPtsU = surf_data["nPtsU"]
        nPtsV = surf_data["nPtsV"]

        # Adding surface to dictionary with tag
        surfaces_data[f"surf_{surfTag}"]=surf_data

        # Updating nodeStartTag
        nodeStartTag = surf_data["nodeEndTag"] + 1

        # Updating surfTag
        surfTag += 1


    if i == range(len(uCuts))[-1]:
        surf_top = cutSurfaces_top[1]
        surf_bottom = cutSurfaces_bottom[1]

        surfList.append(surf_top)
        surfList.append(surf_bottom)

        # Now adding surfaces to the global surfaces dictionary

        surfaces=[surf_top,surf_bottom]

        for surf in surfaces:

           surf_data = {}

           # Flipping control point to u,v format
           surf_data["nodeStartTag"] = nodeStartTag
           surf_data["ctrlpts"] = np.array(
               compatibility.flip_ctrlpts2d(surf.ctrlpts2d[:]))
           surf_data["uKnot"] = surf.knotvector_u
           surf_data["vKnot"] = surf.knotvector_v
           surf_data["p"] = surf.degree_u
           surf_data["q"] = surf.degree_v
           surf_data["nPtsU"] = surf.ctrlpts_size_u
           surf_data["nPtsV"] = surf.ctrlpts_size_v
           surf_data["nodeEndTag"] = surf_data["nodeStartTag"] + \
               surf_data["nPtsU"]*surf_data["nPtsV"] - 1
           surf_data["type"] = "Patch"

           # Adding surface to dictionary with tag
           surfaces_data[f"surf_{surfTag}"]=surf_data

           # Updating nodeStartTag
           nodeStartTag = surf_data["nodeEndTag"] + 1

           # Updating surfTag
           surfTag += 1

    # Updating the remainder surfaces
    topSurf = cutSurfaces_top[1]
    bottomSurf = cutSurfaces_bottom[1]

# VISUALIZATION

surfList_vis = []

for i in range(int(len(surfList)/2)):
    # Translate surfaces for visualization
    surfList_vis.append(operations.translate(
        surfList[2*i], (0, -1e-2, i/40), inplace=False))
    surfList_vis.append(operations.translate(
        surfList[2*i+1], (0, 1e-2, i/40), inplace=False))



container = multi.SurfaceContainer(surfList_vis)
container.sample_size = 20

for surf in container:
    surf.evaluate()


# Visualization configuration
container.vis = VisVTK.VisSurface(ctrlpts=False)

# Render the aspa
container.render()

# surfList now has the surfaces in this order [top_0,bottom_0,top_1,bottom_1, ..... , top_n,bottom_n]

# Bending strip mapping along "u" direction (length)
bendingStripMap_U = []
for i in range(int(len(surfList)/2)):
    bendingStripMap_U.append([2*i, 2*i+1])
    bendingStripMap_U.append([2*i+1, 2*i])

# Bending strip mapping along "v direction (width)
bendingStripMap_V = []
for i in range(int((len(surfList)-2) / 2)):
    bendingStripMap_V.append([2*i, 2*i+2])
    bendingStripMap_V.append([2*i+1, 2*i+3])


# Tag of current bending strip, to recognize the connected patches and its edges
bendingStripTag=0

bendingStrip_dict={}

# Creating bending strips along length

for i in range(len(bendingStripMap_U)):
    bendingStrip_u = bendingStripMap_U[i]

    if (i == 0 or i % 2 == 0):
        edge1 = ["00", "10"]
        edge2 = ["01", "11"]
        retainedSurf = surfList[bendingStrip_u[0]]
        constrainedSurf = surfList[bendingStrip_u[1]]
    else:
        # continue
        edge1 = ["01", "11"]
        edge2 = ["00", "10"]
        constrainedSurf = surfList[bendingStrip_u[0]]
        retainedSurf = surfList[bendingStrip_u[1]]

    # Adding info to bending strips dict local
    bendingStrip_dict_local = {}
    bendingStrip_dict_local["connectedSurfacesTags"] = bendingStrip_u
    bendingStrip_dict_local["connectedEdges"]=[edge1,edge2]

    # Adding info to bending strips dict local
    bendingStrip_dict[f"bendingStrip_{bendingStripTag}"]=bendingStrip_dict_local

    # Updating bendingStripTag
    bendingStripTag+=1

    # Creating first bending strip
    interfacePoints, nodesOnRetainedPatch = edgeGetter(retainedSurf, *edge1)

    nodesOnConstrainedPatch = edgeGetter(constrainedSurf, *edge2)[1]

    bendingStrip = makeBendingStrip(
        nodesOnRetainedPatch, interfacePoints, nodesOnConstrainedPatch, "v")  # From top to bottom



    surfList.append(bendingStrip)

    # Now adding surfaces to the global surfaces dictionary

    surf=bendingStrip

    surf_data = {}

    # Flipping control point to u,v format
    surf_data["nodeStartTag"] = nodeStartTag
    surf_data["ctrlpts"] = np.array(
        compatibility.flip_ctrlpts2d(surf.ctrlpts2d[:]))
    surf_data["uKnot"] = surf.knotvector_u
    surf_data["vKnot"] = surf.knotvector_v
    surf_data["p"] = surf.degree_u
    surf_data["q"] = surf.degree_v
    surf_data["nPtsU"] = surf.ctrlpts_size_u
    surf_data["nPtsV"] = surf.ctrlpts_size_v
    surf_data["nodeEndTag"] = surf_data["nodeStartTag"] + \
        surf_data["nPtsU"]*surf_data["nPtsV"] - 1
    surf_data["type"] = "BendingStrip"

    # Adding surface to dictionary with tag
    surfaces_data[f"surf_{surfTag}"]=surf_data

    # Updating nodeStartTag
    nodeStartTag = surf_data["nodeEndTag"] + 1

    # Updating surfTag
    surfTag += 1


    # Translate surfaces for visualization
    if i == 0 or i % 2 == 0:
        surfList_vis.append(operations.translate(
            bendingStrip, (0, 0, (i/2)/40), inplace=False))
    else:
        surfList_vis.append(operations.translate(
            bendingStrip, (0, 0, (i-1)/2/40), inplace=False))


# Creating bending strips along width
for i in range(len(bendingStripMap_V)):
    bendingStrip_v = bendingStripMap_V[i]

    if (i == 0 or i % 2 == 0):
        edge1 = ["10", "11"]
        edge2 = ["00", "01"]
        retainedSurf = surfList[bendingStrip_v[0]]
        constrainedSurf = surfList[bendingStrip_v[1]]
    else:
        # continue
        edge1 = ["00", "01"]
        edge2 = ["10", "11"]
        constrainedSurf = surfList[bendingStrip_v[0]]
        retainedSurf = surfList[bendingStrip_v[1]]

    # Adding info to bending strips dict local
    bendingStrip_dict_local = {}
    bendingStrip_dict_local["connectedSurfacesTags"] = bendingStrip_v
    bendingStrip_dict_local["connectedEdges"]=[edge1,edge2]

    # Adding info to bending strips dict local
    bendingStrip_dict[f"bendingStrip_{bendingStripTag}"]=bendingStrip_dict_local

    # Updating bendingStripTag
    bendingStripTag+=1

    # Creating first bending strip
    interfacePoints, nodesOnRetainedPatch = edgeGetter(retainedSurf, *edge1)

    nodesOnConstrainedPatch = edgeGetter(constrainedSurf, *edge2)[1]


    bendingStrip = makeBendingStrip(
        nodesOnRetainedPatch, interfacePoints, nodesOnConstrainedPatch, "u")  # From top to bottom



    surfList.append(bendingStrip)

    # Now adding surfaces to the global surfaces dictionary

    surf=bendingStrip

    surf_data = {}

    # Flipping control point to u,v format
    surf_data["nodeStartTag"] = nodeStartTag
    surf_data["ctrlpts"] = np.array(
        compatibility.flip_ctrlpts2d(surf.ctrlpts2d[:]))
    surf_data["uKnot"] = surf.knotvector_u
    surf_data["vKnot"] = surf.knotvector_v
    surf_data["p"] = surf.degree_u
    surf_data["q"] = surf.degree_v
    surf_data["nPtsU"] = surf.ctrlpts_size_u
    surf_data["nPtsV"] = surf.ctrlpts_size_v
    surf_data["nodeEndTag"] = surf_data["nodeStartTag"] + \
        surf_data["nPtsU"]*surf_data["nPtsV"] - 1
    surf_data["type"] = "BendingStrip"

    # Adding surface to dictionary with tag
    surfaces_data[f"surf_{surfTag}"]=surf_data

    # Updating nodeStartTag
    nodeStartTag = surf_data["nodeEndTag"] + 1

    # Updating surfTag
    surfTag += 1

    # Translate surfaces for visualization
    if i == 0 or i % 2 == 0:
        surfList_vis.append(operations.translate(
            bendingStrip, (0, -2e-2, (i/2)/40), inplace=False))
    else:
        surfList_vis.append(operations.translate(
            bendingStrip, (0, 2e-2, (i-1)/2/40), inplace=False))


# TO-DO: define bending strips according to this order

# surfList.append(topSurf)
# surfList.append(bottomSurf)

container = multi.SurfaceContainer(surfList_vis)
container.sample_size = 20

for surf in container:
    surf.evaluate()


# Visualization configuration
container.vis = VisVTK.VisSurface(ctrlpts=False)

# Render the aspa
container.render()


# Original surface
container = multi.SurfaceContainer(surfList)
container.sample_size = 20

for surf in container:
    surf.evaluate()


# Visualization configuration
container.vis = VisVTK.VisSurface(ctrlpts=True)

# Render the aspa
evalcolor = ["green", "red", "blue", "black"]
# container.render(evalcolor=evalcolor)
container.render()












save_dict(surfaces_data, 'surfaces_data.pkl')
save_dict(bendingStrip_dict, 'bendingStrip_dict.pkl')












# print(surf01_data_loaded)


exit()


# Creating first bending strip
interfacePoints, nodesOnRightPatch = edgeGetter(topSurf, "01", "11")

nodesOnLeftPatch = edgeGetter(bottomSurf, "00", "10")[1]

bendingStrip_1 = makeBendingStrip(
    nodesOnRightPatch, interfacePoints, nodesOnLeftPatch)  # From top to bottom

# Creating second bending strip
interfacePoints, nodesOnRightPatch = edgeGetter(topSurf, "00", "10")

nodesOnLeftPatch = edgeGetter(bottomSurf, "01", "11")[1]

bendingStrip_2 = makeBendingStrip(
    nodesOnRightPatch, interfacePoints, nodesOnLeftPatch)  # From top to bottom

# exit()
# bendingStrip_1 = surfFromFile("coarse.3.dat")
# bendingStrip_2 = surfFromFile("coarse.4.dat")

surfList.append(topSurf)
surfList.append(bottomSurf)
surfList.append(bendingStrip_1)
surfList.append(bendingStrip_2)


container = multi.SurfaceContainer(surfList)
container.sample_size = 20

for surf in container:
    surf.evaluate()


# Visualization configuration
container.vis = VisVTK.VisSurface(ctrlpts=True, legend=False)

# Render the aspa
evalcolor = ["green", "red", "blue", "blue"]
container.render(evalcolor=evalcolor)


ops.wipe()
ops.model('basic', '-ndm', 3, '-ndf', 3)
# E_prom=(30300e6+58300e6)*0.5
# E_prom=58300e6


# Real Data
E1 = 58300e6  # Young's modulus N/m^2 Pa
E2 = 58300e6
nu12 = 0.055  # Poisson's ratio
nu21 = 0.055  # Poisson's ratio
G12 = 4118.3e6  # Shear modulus
rho = 2900  # *9.807 # kg/m^3
mm = 1.0 / 1000
t = 0.005


# # Testing Data
# E1 = 2.1e11  # Young's modulus N/m^2 Pa
# E2 = 2.1e11
# nu12 = 0.3 # Poisson's ratio
# nu21 = 0.3  # Poisson's ratio
# G12 = E1/(2*(1+nu12))  # Shear modulus
# rho = 2900  # *9.807 # kg/m^3
# mm = 1.0 / 1000
# t = 0.05


# tagNDmat1 = 1
# ops.nDMaterial("ElasticIsotropic", tagNDmat1, E1, nu12, rho)

# # nDMaterial PlateFiber $nDtag_platefiber $nDtag_elastic
# tagPlaneStress1 = 2
# ops.nDMaterial("PlaneStress", tagPlaneStress1, tagNDmat1)

tagPlaneStress1 = 1
ops.nDMaterial("ElasticOrthotropicPlaneStress",
               tagPlaneStress1, E1, E2, nu12, nu21, G12, rho)
tagPlaneStress2 = 2
ops.nDMaterial("ElasticOrthotropicPlaneStress",
               tagPlaneStress2, E1, E2, nu12, nu21, G12, rho)

tagPlaneStress3 = 3
ops.nDMaterial("ElasticOrthotropicPlaneStress", tagPlaneStress3,
               0*E1 * 1e5, E1*1e5, 0, 0, 0, rho)
tagPlaneStress4 = 4
ops.nDMaterial("ElasticOrthotropicPlaneStress", tagPlaneStress4,
               0*E1 * 1e5, E1*1e5, 0, 0, 0, rho)

deg2rad = pi / 180

matTags = [1]
thickness = [t]
θ = [0 * deg2rad]


g = -9.807
gFact = [0, 0*g, 0]


Nlayers = len(θ)

names = ["topPatch", "bottomPatch", "bendingStrip_01", "bendingStrip_02"]
shellType = 'KLShell'

patchTag = 1
nodeStartTag = 1

nodesMap = []
patchTags = []

for i in range(len(container)):
    patchTags.append(patchTag)
    surf = container[i]
    name = names[i]

    # Identifying bending strip
    if name in ["bendingStrip_01", "bendingStrip_02"]:
        shellType = "KLShell_BendingStrip"

    # Flipping control point to u,v format
    controlPts = surf.ctrlpts2d[:]
    controlPts = np.array(compatibility.flip_ctrlpts2d(controlPts))

    for dim in controlPts:  # Unweighting control pts
        for point in dim:
            point[0:3] /= point[3]

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
    nodesMap.append(
        np.arange(nodeStartTag, ops.getNodeTags()[-1] + 1).tolist())

    # Update patchTag, nodeStartTag and materialTag
    lastElTag = ops.getEleTags()[-1]
    lastNodeTag = ops.getNodeTags()[-1]
    matTags[0] += 1
    patchTag = lastElTag + 1
    nodeStartTag = lastNodeTag + 1


# Getting TopSurf info
print("\nGetting topSurf info\n")
nPoints_topSurf = topSurf.ctrlpts_size_u * topSurf.ctrlpts_size_v

nodesOn_0010_topSurf = np.arange(1, topSurf.ctrlpts_size_u + 1, 1)
nextTo_0010_topSurf = nodesOn_0010_topSurf + topSurf.ctrlpts_size_u

nodesOn_0001_topSurf = np.arange(
    1, nPoints_topSurf + 1, topSurf.ctrlpts_size_u)
nextTo_0001_topSurf = nodesOn_0001_topSurf + 1

nodesOn_0111_topSurf = np.arange(
    nPoints_topSurf - topSurf.ctrlpts_size_u + 1, nPoints_topSurf + 1, 1)
nextTo_0111_topSurf = nodesOn_0111_topSurf - topSurf.ctrlpts_size_u

nodesOn_1011_topSurf = np.arange(
    topSurf.ctrlpts_size_u, nPoints_topSurf + 1, topSurf.ctrlpts_size_u)
nextTo_1011_topSurf = nodesOn_1011_topSurf - 1


# Getting BottomSurf info
print("\nGetting bottomSurf info\n")
nPoints_bottomSurf = bottomSurf.ctrlpts_size_u * bottomSurf.ctrlpts_size_v

nodesOn_0010_bottomSurf = np.arange(
    1, bottomSurf.ctrlpts_size_u + 1, 1) + nPoints_topSurf
nextTo_0010_bottomSurf = nodesOn_0010_bottomSurf + bottomSurf.ctrlpts_size_u

nodesOn_0001_bottomSurf = np.arange(
    1, nPoints_bottomSurf + 1, bottomSurf.ctrlpts_size_u) + nPoints_topSurf
nextTo_0001_bottomSurf = nodesOn_0001_bottomSurf + 1

nodesOn_0111_bottomSurf = np.arange(
    nPoints_bottomSurf - bottomSurf.ctrlpts_size_u + 1, nPoints_bottomSurf + 1, 1) + nPoints_topSurf
nextTo_0111_bottomSurf = nodesOn_0111_bottomSurf - bottomSurf.ctrlpts_size_u

nodesOn_1011_bottomSurf = np.arange(
    bottomSurf.ctrlpts_size_u, nPoints_bottomSurf + 1, bottomSurf.ctrlpts_size_u) + nPoints_topSurf
nextTo_1011_bottomSurf = nodesOn_1011_bottomSurf - 1


print("Creating constraints\n")
# Creating constraints

# # Restraining all nodes on X
# for n in ops.getNodeTags():
#     n=int(n)
#     ops.fix(n,1,0,0)


# fixedNodes=np.concatenate([nodesOn_0001_bottomSurf,nodesOn_0001_topSurf])
# for n in fixedNodes:
#     n=int(n)
#     ops.fix(n,1,1,1)
ops.fixZ(-0.229, 1, 1, 1, '-tol', 1e-5)  # First row
# ops.fixZ(-0.119, 1, 1, 1)  # Second row


print("Equal Dofing\n")
# Equal dofing bending strips and patch overlaps


# First patch overlaps


# bottomSurf,"00","10"
# topSurf,"01","11"

for i in range(len(nodesOn_0111_topSurf)):
    retainedNode = int(nodesOn_0111_topSurf[i])
    constrainedNode = int(nodesOn_0010_bottomSurf[i])
    ops.equalDOF(retainedNode, constrainedNode, 1, 2, 3)


# bottomSurf,"01","11"
# topSurf,"00","10"

for i in range(len(nodesOn_0010_topSurf)):
    retainedNode = int(nodesOn_0010_topSurf[i])
    constrainedNode = int(nodesOn_0111_bottomSurf[i])
    ops.equalDOF(retainedNode, constrainedNode, 1, 2, 3)


# EqualDofing bendingStrip_01

# Retaining Nodes
firstRow_r = nextTo_0111_topSurf
secondRow_r = nodesOn_0111_topSurf
thirdRow_r = nextTo_0010_bottomSurf

# # Constrained nodes
nodesOnBendingStrip_01 = nodesMap[2]

firstRow_c = []
secondRow_c = []
thirdRow_c = []


firstRow_c = (nodesOnBendingStrip_01[:len(firstRow_r)])
secondRow_c = (nodesOnBendingStrip_01[len(firstRow_r):2*len(firstRow_r)])
thirdRow_c = (nodesOnBendingStrip_01[2*len(firstRow_r):])


for i in range(len(firstRow_c)):
    retainedNode = int(firstRow_r[i])
    constrainedNode = int(firstRow_c[i])
    ops.equalDOF(retainedNode, constrainedNode, 1, 2, 3)


for i in range(len(secondRow_c)):
    retainedNode = int(secondRow_r[i])
    constrainedNode = int(secondRow_c[i])
    ops.equalDOF(retainedNode, constrainedNode, 1, 2, 3)


for i in range(len(thirdRow_c)):
    retainedNode = int(thirdRow_r[i])
    constrainedNode = int(thirdRow_c[i])
    ops.equalDOF(retainedNode, constrainedNode, 1, 2, 3)


# EqualDofing bendingStrip_02

# Retaining Nodes
firstRow_r = nextTo_0010_topSurf
secondRow_r = nodesOn_0010_topSurf
thirdRow_r = nextTo_0111_bottomSurf

# # Constrained nodes
nodesOnBendingStrip_02 = nodesMap[3]

firstRow_c = (nodesOnBendingStrip_02[:len(firstRow_r)])
secondRow_c = (nodesOnBendingStrip_02[len(firstRow_r):2*len(firstRow_r)])
thirdRow_c = (nodesOnBendingStrip_02[2*len(firstRow_r):])


for i in range(len(firstRow_c)):
    retainedNode = int(firstRow_r[i])
    constrainedNode = int(firstRow_c[i])
    ops.equalDOF(retainedNode, constrainedNode, 1, 2, 3)


for i in range(len(secondRow_c)):
    retainedNode = int(secondRow_r[i])
    constrainedNode = int(secondRow_c[i])
    ops.equalDOF(retainedNode, constrainedNode, 1, 2, 3)


for i in range(len(thirdRow_c)):
    retainedNode = int(thirdRow_r[i])
    constrainedNode = int(thirdRow_c[i])
    ops.equalDOF(retainedNode, constrainedNode, 1, 2, 3)


print('Finished equalDofing')
# Finished equalDofing


print("\n\n\nPRINTING DOMAIN-----------------------")
ops.printModel()
print("\n\n\nDONE PRINTING DOMAIN-----------------------")


# # ------------------------------
# # Start of analysis generation
# # ------------------------------

# create TimeSeries
ops.timeSeries("Linear", 1)

# create a plain load pattern
ops.pattern("Plain", 1, 1)


print("Loading nodes")
P = 5
# for n in nodesMap[0]:
loadedNodes = []
# loadedNodes=np.concatenate([nodesOn_1011_topSurf,nodesOn_1011_bottomSurf,nextTo_1011_topSurf,nextTo_1011_bottomSurf])
# loadedNodes=np.concatenate([nodesOn_0010_topSurf,nodesOn_0010_bottomSurf, nodesOn_0111_topSurf,nodesOn_0111_bottomSurf])
for n in loadedNodes:
    n = int(n)
    print("n: ", n)
    print("ops.nodeCoord(n): ", ops.nodeCoord(n))
    ops.load(n, 0, P, 0)


angle = 0*deg2rad
weight = [0.0, np.cos(angle)*g, np.sin(angle)*g]
ops.eleLoad("-ele", patchTags[0], "-type", "-SelfWeight", *weight)
ops.eleLoad("-ele", patchTags[1], "-type", "-SelfWeight", *weight)


print("Finished loading nodes")


print("Starting analysis")

# Create test
ops.test("NormDispIncr", 1.0e-4, 30, 1)


# create SOE
# ops.system("FullGeneral")
ops.system("UmfPack")

# create DOF number
# ops.numberer("Plain")
ops.numberer("RCM")

# create constraint handler
ops.constraints("Plain")
# ops.constraints("Lagrange")
# ops.constraints("Transformation")
# ops.constraints("Penalty",1e10,1e10)


ops.algorithm("Linear")
# ops.algorithm("Newton")
# ops.algorithm("SecantNewton")
# ops.algorithm("NewtonLineSearch", 'type', 'Bisection')
# ops.algorithm("NewtonLineSearch")
# ops.algorithm("ModifiedNewton")
# ops.algorithm("KrylovNewton")
# ops.algorithm("BFGS")
# ops.algorithm("Broyden")

# create integrator

nSteps = 1
ops.integrator("LoadControl", 1.0 / nSteps)


# create analysis object
ops.analysis("Static")


for j in range(nSteps):
    print("=================================")
    print(f"Load step {j}")
    print("=================================")
    result = ops.analyze(1)
    if result != 0:
        break
        exit(-1)
    else:
        # Adding deformation to control points
        fdef = 1e2

        for i in range(len(container)):
            surf = container[i]
            nodes = nodesMap[i]
            controlPts = np.array(compatibility.flip_ctrlpts2d(surf.ctrlpts2d[:])).reshape(
                [surf.ctrlpts_size_v * surf.ctrlpts_size_u, 4]).tolist()
            for n in nodes:
                # Get node position in list
                indexN = nodes.index(n)
                point = controlPts[indexN]

                # Add deformation scaled by fdef
                weight = point[3]
                point[0] = (ops.nodeCoord(n, 1) + fdef *
                            ops.nodeDisp(n, 1)) * weight
                point[1] = (ops.nodeCoord(n, 2) + fdef *
                            ops.nodeDisp(n, 2)) * weight
                point[2] = (ops.nodeCoord(n, 3) + fdef *
                            ops.nodeDisp(n, 3)) * weight

            nPoints = surf.ctrlpts_size_u * surf.ctrlpts_size_v
            shape = np.array(compatibility.flip_ctrlpts2d(
                surf.ctrlpts2d[:])).shape
            controlPts = np.array(controlPts).reshape(shape)
            controlPts = np.array(compatibility.flip_ctrlpts2d(controlPts))

            surf.set_ctrlpts(controlPts.reshape(nPoints, 4).tolist(),
                             surf.ctrlpts_size_u, surf.ctrlpts_size_v)

        # Visualizing deformed surface
        container.sample_size = 40
        for surf in container:
            surf.evaluate()

        # Visualization configuration
        container.vis = VisVTK.VisSurface(ctrlpts=False, legend=False)

        # Render the aspa
        container.render(evalcolor=evalcolor)

        for i in range(len(container)):
            surf = container[i]
            nodes = nodesMap[i]
            controlPts = np.array(compatibility.flip_ctrlpts2d(surf.ctrlpts2d[:])).reshape(
                [surf.ctrlpts_size_v * surf.ctrlpts_size_u, 4]).tolist()
            for n in nodes:
                # Get node position in list
                indexN = nodes.index(n)
                point = controlPts[indexN]

                # Add deformation scaled by fdef
                weight = point[3]
                point[0] = (ops.nodeCoord(n, 1)) * weight
                point[1] = (ops.nodeCoord(n, 2)) * weight
                point[2] = (ops.nodeCoord(n, 3)) * weight

            nPoints = surf.ctrlpts_size_u * surf.ctrlpts_size_v
            shape = np.array(compatibility.flip_ctrlpts2d(
                surf.ctrlpts2d[:])).shape
            controlPts = np.array(controlPts).reshape(shape)
            controlPts = np.array(compatibility.flip_ctrlpts2d(controlPts))

            surf.set_ctrlpts(controlPts.reshape(nPoints, 4).tolist(),
                             surf.ctrlpts_size_u, surf.ctrlpts_size_v)

        print("\nNext load step\n")


# exit()

# Starting manual eigen
fixedNodes = np.array(nodesOn_0001_topSurf)-1
noCtrPts = len(ops.getNodeTags())
noDofs = 3*noCtrPts
freeNodes = list(set(np.arange(noCtrPts))-set(fixedNodes))


# % essential boundary conditions (fixed along edge)
udofs = 3*(fixedNodes+1)-2 - 1
vdofs = 3*(fixedNodes+1)-1 - 1
wdofs = 3*(fixedNodes+1)-0 - 1


uFixed = np.zeros((len(fixedNodes)))
vFixed = np.zeros((len(fixedNodes)))
wFixed = np.zeros((len(fixedNodes)))

oldDofs = range(noDofs)
freeDofs = (list(set(oldDofs) - set(udofs)-set(vdofs)-set(wdofs)))


# Getting matrices from GimmeMCK
ops.wipeAnalysis()
ops.system('FullGeneral')
ops.analysis('Transient')

# Mass
ops.integrator('GimmeMCK', 1.0, 0.0, 0.0)
ops.analyze(1, 0.0)

# Number of equations in the model
N = ops.systemSize()  # Has to be done after analyze

M = ops.printA('-ret')  # Or use ops.printA('-file','M.out')
M = np.array(M)  # Convert the list to an array
M.shape = (N, N)  # Make the array an NxN matrix


# Stiffness
ops.integrator('GimmeMCK', 0.0, 0.0, 1.0)
ops.analyze(1, 0.0)
K = ops.printA('-ret')
K = np.array(K)
K.shape = (N, N)


# Kff = K[np.ix_(freeDofs, freeDofs)]
# Mff = M[np.ix_(freeDofs, freeDofs)]
# W, phi_full = eigsh(Kff, M=Mff, k=10, sigma=0, maxiter=1000000000)

W, phi_full = eigsh(K, M=M, k=15, sigma=0, maxiter=1000000000)

W = np.sqrt(W)
order = np.argsort(W)
W = W[order].real

print("W: ", W)
print("W/(2*np.pi): ", W/(2*np.pi))
# [ 129.92415927  315.65654987  387.76783551  504.37360884  566.74681745
#   708.16318993  739.76512482  924.18804355 1021.92342721 1230.60220003]


# exit()


# Eigenvalues

print("Starting Eigen")

# ops.wipeAnalysis()

# # create SOE
# ops.system("FullGeneral")

# # create DOF number
# ops.numberer("Plain")

# # create constraint handler
# ops.constraints("Plain")

# ops.system("BandSPD")
# ops.integrator("Newmark", 0.5, 0.25)
# ops.system("FullGeneral")

nodes = ops.getNodeTags()

Nnodes = len(nodes)
Neigenvalues = 10  # arpack can only compute N-1 eigvals


w2s = ops.eigen(Neigenvalues)
# w2s = ops.eigen('-solver','-fullGenLapack',Neigenvalues)
# w2s = ops.eigen('-standard','-symmBandLapack',Neigenvalues)


order = np.argsort(w2s)
w2s = np.array(w2s, dtype=np.float64)[order]
w = np.sqrt(w2s)
print("w: ", w)
# exit()
# exit()


for i, w2 in enumerate(w2s):
    w = sqrt(abs(w2))
    f = w / 2 / pi
    T = 1 / f
    print(f"{i} {w2} {w} {f} {T} ")


phi = np.zeros((3 * Nnodes, Neigenvalues))

for i, n in enumerate(nodes):
    for j in range(Neigenvalues):
        phi[3 * i:3 * i + 3, j] = ops.nodeEigenvector(n, j + 1)

print(f"ϕ = {phi}")

w = np.sqrt(abs(w2s))
print("w/2/pi: ", w/2/pi)
for j in range(Neigenvalues):
    w = sqrt(abs(w2s[j]))
    f = w / 2 / pi
    print("=================================")
    print(f"Eigenvalue {j}")
    print(f"f = {f} Hz")
    print("=================================")

    # Adding deformation to control points
    fdef = 1e0
    for i in range(len(container)):
        surf = container[i]
        nodes = nodesMap[i]
        controlPts = np.array(compatibility.flip_ctrlpts2d(surf.ctrlpts2d[:])).reshape(
            [surf.ctrlpts_size_v * surf.ctrlpts_size_u, 4]).tolist()
        for n in nodes:
            # Get node position in list
            indexN = nodes.index(n)
            point = controlPts[indexN]

            eigenVector = phi[:, order[j]]

            # Add deformation scaled by fdef
            weight = point[3]
            # point[0] = (ops.nodeCoord(n)[0] + fdef * eigenVector[3*(n-1)] )* weight
            # point[1] = (ops.nodeCoord(n)[1] + fdef * eigenVector[3*(n-1)+1] )* weight
            # point[2] = (ops.nodeCoord(n)[2] + fdef * eigenVector[3*(n-1)+2] )* weight
            # print("n: ", n)
            # print("j: ", j)
            # print("ops.nodeEigenvector(n,j): ", ops.nodeEigenvector(n,j+1))

            # From portwood digital
            point[0] = (ops.nodeCoord(n)[0] + fdef *
                        ops.nodeEigenvector(n, j+1, 1)) * weight
            point[1] = (ops.nodeCoord(n)[1] + fdef *
                        ops.nodeEigenvector(n, j+1, 2)) * weight
            point[2] = (ops.nodeCoord(n)[2] + fdef *
                        ops.nodeEigenvector(n, j+1, 3)) * weight

        nPoints = surf.ctrlpts_size_u * surf.ctrlpts_size_v
        shape = np.array(compatibility.flip_ctrlpts2d(surf.ctrlpts2d[:])).shape
        controlPts = np.array(controlPts).reshape(shape)
        controlPts = np.array(compatibility.flip_ctrlpts2d(controlPts))

        surf.set_ctrlpts(controlPts.reshape(nPoints, 4).tolist(),
                         surf.ctrlpts_size_u, surf.ctrlpts_size_v)

    # Visualizing deformed surface
    container.sample_size = 40
    for surf in container:
        surf.evaluate()

    # Render the aspa
    container.render(evalcolor=evalcolor)

    for i in range(len(container)):
        surf = container[i]
        nodes = nodesMap[i]
        controlPts = np.array(compatibility.flip_ctrlpts2d(surf.ctrlpts2d[:])).reshape(
            [surf.ctrlpts_size_v * surf.ctrlpts_size_u, 4]).tolist()
        for n in nodes:
            # Get node position in list
            indexN = nodes.index(n)
            point = controlPts[indexN]

            # Add deformation scaled by fdef
            weight = point[3]
            point[0] = (ops.nodeCoord(n)[0]) * weight
            point[1] = (ops.nodeCoord(n)[1]) * weight
            point[2] = (ops.nodeCoord(n)[2]) * weight

        nPoints = surf.ctrlpts_size_u * surf.ctrlpts_size_v
        shape = np.array(compatibility.flip_ctrlpts2d(surf.ctrlpts2d[:])).shape
        controlPts = np.array(controlPts).reshape(shape)
        controlPts = np.array(compatibility.flip_ctrlpts2d(controlPts))

        surf.set_ctrlpts(controlPts.reshape(nPoints, 4).tolist(),
                         surf.ctrlpts_size_u, surf.ctrlpts_size_v)

    print("\nNext Eigenvalue\n")
