# Function to handle edges. Receives surface and returns desired edge and nodes next to edge


import numpy as np

# Geomgl utilities for visualization and surface manipulation
from geomdl import NURBS, compatibility, operations, multi, knotvector


def edgeGetter(surf, pointA, pointB):
    controlPts_vu = surf.ctrlpts2d[::]  # These are given in v,u
    controlPts_uv = compatibility.flip_ctrlpts2d(controlPts_vu)  # These are given in v,u
    # print("controlPts_vu: ", controlPts_vu)
    # print("controlPts_uv: ", controlPts_uv)

    if pointA == "10" and pointB == "11":
        edgePoints = controlPts_vu[-1]
        nextToEdgePoints = controlPts_vu[-2]

    elif pointA == "00" and pointB == "10":
        edgePoints = controlPts_uv[0]
        nextToEdgePoints = controlPts_uv[1]

    elif pointA == "01" and pointB == "11":
        edgePoints = controlPts_uv[-1]
        nextToEdgePoints = controlPts_uv[-2]

    elif pointA == "00" and pointB == "01":
        edgePoints = controlPts_vu[0]
        nextToEdgePoints = controlPts_vu[1]

    else:
    	print("Unexistant edge")
    	edgePoints=[]
    	nextToEdgePoints=[]

    return edgePoints,nextToEdgePoints


# def makeBendingStrip(pointsL,interface,pointsR):
# 	# Receives u,v, returns bending strip surface

# 	if len(pointsL)!=len(pointsR)!=len(interface):
# 		print("Non matching dimensions")
# 		return []

# 	sizeU=3
# 	sizeV=len(pointsL)

# 	controlPts=np.zeros([sizeU*sizeV,4])

# 	k=0
# 	for i in range(sizeV):
# 		controlPts[k,:]=pointsL[i]
# 		k+=1
# 	for i in range(sizeV):
# 		controlPts[k,:]=interface[i]
# 		k+=1
# 	for i in range(sizeV):
# 		controlPts[k,:]=pointsR[i]
# 		k+=1

# 	# Create surface 
# 	bendingStrip = NURBS.Surface()

# 	# Set surface degrees
# 	bendingStrip.degree_u = 2
# 	bendingStrip.degree_v = 1


# 	# for point in controlPts:
# 	# 	# point[0:3] /= point[3]
# 	# 	point[0:3] *= point[3]
# 		# point[3] = 1

# 	# Setting control points for surface
# 	bendingStrip.set_ctrlpts(controlPts.tolist(), sizeU, sizeV)

# 	# Set knot vectors
# 	bendingStrip.knotvector_u = knotvector.generate(bendingStrip.degree_u, bendingStrip.ctrlpts_size_u)
# 	bendingStrip.knotvector_v = knotvector.generate(bendingStrip.degree_v, bendingStrip.ctrlpts_size_v)

# 	return bendingStrip

def makeBendingStrip(pointsL,interface,pointsR,direction):
	# Receives u,v, returns bending strip surface

	if len(pointsL)!=len(pointsR)!=len(interface):
		print("Non matching dimensions")
		return []

	if direction == 'u':
		sizeU=3
		sizeV=len(pointsL)
	else:
		sizeU=len(pointsL)
		sizeV=3

	controlPts=np.zeros([sizeU*sizeV,4])

	if direction == "u":
		k=0
		for i in range(sizeV):
			controlPts[i,:]=pointsL[k]
			k+=1
		k=0
		for i in range(sizeV,2*sizeV):
			controlPts[i,:]=interface[k]
			k+=1
		k=0
		for i in range(2*sizeV,3*sizeV):
			controlPts[i,:]=pointsR[k]
			k+=1

	else:
		for i in range(sizeU):
			controlPts[3*i,:]=pointsL[i]
			controlPts[3*i+1,:]=interface[i]
			controlPts[3*i+2,:]=pointsR[i]



	# Create surface 
	bendingStrip = NURBS.Surface()

	# Set surface degrees
	if direction == "v":
		bendingStrip.degree_u = 1
		bendingStrip.degree_v = 2
	else:
		bendingStrip.degree_u = 2
		bendingStrip.degree_v = 1



	# for point in controlPts:
	# 	# point[0:3] /= point[3]
	# 	point[0:3] *= point[3]
		# point[3] = 1

	# Setting control points for surface
	bendingStrip.set_ctrlpts(controlPts.tolist(), sizeU, sizeV)

	# Set knot vectors
	bendingStrip.knotvector_u = knotvector.generate(bendingStrip.degree_u, bendingStrip.ctrlpts_size_u)
	bendingStrip.knotvector_v = knotvector.generate(bendingStrip.degree_v, bendingStrip.ctrlpts_size_v)

	return bendingStrip
