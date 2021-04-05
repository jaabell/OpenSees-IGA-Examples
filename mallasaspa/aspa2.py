import os
from geomdl import NURBS, multi
import numpy as np
from surfVisualize import surfVisualize
from geomdl import operations, compatibility, exchange, exchange_vtk



# Fix file path
os.chdir(os.path.dirname(os.path.realpath(__file__)))

from fileGetter import surfFromFile, generateKnotVector

surfList=[]
topSurf=surfFromFile("coarse.1.dat",deg_elevate=True)
# exit()
bottomSurf=surfFromFile("coarse.2.dat",deg_elevate=True)
bendingStrip_1=surfFromFile("coarse.3.dat")
bendingStrip_2=surfFromFile("coarse.4.dat")


surfList.append(topSurf)
surfList.append(bottomSurf)
surfList.append(bendingStrip_1)
surfList.append(bendingStrip_2)

aspa=multi.SurfaceContainer(surfList)
aspa.sample_size=25

for surf in surfList:
	operations.refine_knotvector(surf, [0, 0])
	surf.knotvector_u=generateKnotVector(surf.degree_u, surf.ctrlpts_size_u)
	surf.knotvector_v=generateKnotVector(surf.degree_v, surf.ctrlpts_size_v)

for surf in aspa:
	surf.evaluate()


exchange.export_stl(topSurf, "top.stl",vertex_spacing=1)
# Export evaluated points as a .vtk file
exchange_vtk.export_polydata(topSurf, "top.vtk", tessellate=True)

exchange.export_stl(bottomSurf, "bottom.stl",vertex_spacing=1)
exchange.export_stl(bendingStrip_1, "bendingStrip_1.stl")
exchange.export_stl(bendingStrip_2, "bendingStrip_2.stl")

# print("secondRing: ", secondRing)




# for i in range(25):
# 	print("ctrlpts[i]: ", ctrlpts[i])
# bottomSurf_ctrlPts = compatibility.flip_ctrlpts2d(bottomSurf.ctrlpts2d)


# print("bottomSurf_ctrlPts: ", bottomSurf_ctrlPts)


# print("bottomSurf_ctrlPts[:][0]: ", bottomSurf_ctrlPts[:][1])
# for i in range(24):
# 	bottomSurf_ctrlPts[:][0][i][3]=10

# print("bottomSurf_ctrlPts: ", bottomSurf_ctrlPts)



# Refine knot vectors and update geometry






from geomdl.visualization import VisMPL

# Visualization configuration
aspa.vis=VisMPL.VisSurface(ctrlpts=False, legend=False, figure_size=[940, 940])

# Render the ducky
evalcolor=["cyan",'cyan','yellow','yellow']
aspa.render(evalcolor=evalcolor)


# Surface visualization with geomdl functions and multiprocessing for holding

# surfVisualize(topSurf,visualizer="VTK",sample_size=200)
# exit()
# surfVisualize(surfaces,visualizer="MPL",sample_size=50)
