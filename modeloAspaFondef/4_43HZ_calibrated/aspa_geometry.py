
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


# Input parameters
refine_U = True
uCuts = np.array([0.011, 0.251, 0.696, 1.141, 1.586])  # Places to define sections
knotInsert_first = 0.05
knotInsert_last = 1 - knotInsert_first
knot_insert = True

# For each surface, I need it's control points, order, knot vectors, startNodeTag, endNodeTag, type
surfaces_data = {}

surfList = []
topSurf_original = surfFromFile("coarse.1.dat",refine_U)
bottomSurf_original = surfFromFile("coarse.2.dat",refine_U)
# 

# Interpolating length, in order to know where to split the surfaces

# Getting the bounding box, average between two surfaces
boundingBox = (np.array(topSurf_original.bbox) +
               np.array(bottomSurf_original.bbox))/2
zLimits = boundingBox[:, 2]  # Getting the limits for parametrization
zLength = sum(abs(zLimits))

parametricInt = interpolate.interp1d(zLimits, [0, 1])

# Parametric coordinates (along "u") that define the sections
parametricCuts = parametricInt(uCuts)


# Adding a knot where the furthest accelerometer is : 1.991 along length
# u_accelerometer = float(parametricInt(1.991))
# operations.insert_knot(topSurf_original,[u_accelerometer,None],[1,0])
# operations.insert_knot(bottomSurf_original,[u_accelerometer,None],[1,0])

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
  
    surf_top = cutSurfaces_top[0]
    surf_bottom = cutSurfaces_bottom[0]


    if knot_insert:
        # Adding knots so the bending strips are not to far inside the patch
        operations.insert_knot(surf_top,[knotInsert_first,None],[1,0])
        operations.insert_knot(surf_bottom,[knotInsert_first,None],[1,0])

        operations.insert_knot(surf_top,[knotInsert_last,None],[1,0])
        operations.insert_knot(surf_bottom,[knotInsert_last,None],[1,0])

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
        surf_data["zone"] = i


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

        if knot_insert:
            # Adding knots so the bending strips are not to far inside the patch
            operations.insert_knot(surf_top,[knotInsert_first,None],[1,0])
            operations.insert_knot(surf_bottom,[knotInsert_first,None],[1,0])

            operations.insert_knot(surf_top,[knotInsert_last,None],[1,0])
            operations.insert_knot(surf_bottom,[knotInsert_last,None],[1,0])


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
           surf_data["zone"] = i + 1

        

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
bendingStripTag=surfTag

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

    # Adding info to bending strips dict local
    bendingStrip_dict_local["surfTag"]=f"surf_{surfTag}"
    bendingStrip_dict[f"bendingStrip_{bendingStripTag}"]=bendingStrip_dict_local

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
        edge1 = ["10", "11"]
        edge2 = ["00", "01"]
        constrainedSurf = surfList[bendingStrip_v[1]]
        retainedSurf = surfList[bendingStrip_v[0]]

    # Adding info to bending strips dict local
    bendingStrip_dict_local = {}
    bendingStrip_dict_local["connectedSurfacesTags"] = bendingStrip_v
    bendingStrip_dict_local["connectedEdges"]=[edge1,edge2]




    # Adding info to bending strips dict global
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

