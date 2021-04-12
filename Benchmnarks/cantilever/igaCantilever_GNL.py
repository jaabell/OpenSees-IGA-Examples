#  IGA CANTILEVER PLATE with bending strip


import numpy as np
import opensees as ops
from math import *

# Geomgl utilities for visualization and surface manipulation
from geomdl import NURBS, compatibility, operations, multi, knotvector
from surfVisualize import *

La = 5.0  	#
Lb = 1.0  	#
mm = 1.0 / 1000.  # m

ops.wipe()
ops.model('basic', '-ndm', 3, '-ndf', 3)

# These are given in v,u
controlPts = np.array([
  [0.00000  , 0.0 , 0.0 , 1.0] ,
  [0.00000  , 1.0 , 0.0 , 1.0] ,
  [0.31250  , 0.0 , 0.0 , 1.0] ,
  [0.31250  , 1.0 , 0.0 , 1.0] ,
  [0.93750  , 0.0 , 0.0 , 1.0] ,
  [0.93750  , 1.0 , 0.0 , 1.0] ,
  [1.87500  , 0.0 , 0.0 , 1.0] ,
  [1.87500  , 1.0 , 0.0 , 1.0] ,
  [3.12500  , 0.0 , 0.0 , 1.0] ,
  [3.12500  , 1.0 , 0.0 , 1.0] ,
  [4.37500  , 0.0 , 0.0 , 1.0] ,
  [4.37500  , 1.0 , 0.0 , 1.0] ,
  [5.62500  , 0.0 , 0.0 , 1.0] ,
  [5.62500  , 1.0 , 0.0 , 1.0] ,
  [6.87500  , 0.0 , 0.0 , 1.0] ,
  [6.87500  , 1.0 , 0.0 , 1.0] ,
  [8.12500  , 0.0 , 0.0 , 1.0] ,
  [8.12500  , 1.0 , 0.0 , 1.0] ,
  [9.06250  , 0.0 , 0.0 , 1.0] ,
  [9.06250  , 1.0 , 0.0 , 1.0] ,
  [9.68750  , 0.0 , 0.0 , 1.0] ,
  [9.68750  , 1.0 , 0.0 , 1.0] ,
  [10.00000 , 0.0 , 0.0 , 1.0] ,
  [10.00000 , 1.0 , 0.0 , 1.0]
])




# Create a NURBS surface instance
surf = NURBS.Surface()

# Set surface degrees
surf.degree_u = 4
surf.degree_v = 1

# Setting control points for surface
surf.set_ctrlpts(controlPts.tolist(), 12, 2)

# Set knot vectors
surf.knotvector_u = knotvector.generate(surf.degree_u, surf.ctrlpts_size_u)
surf.knotvector_v = knotvector.generate(surf.degree_v, surf.ctrlpts_size_v)


# nDMaterial ElasticIsotropic $nDtag_elastic $elasticidad_probeta
# $poisson_probeta
E1 = 1.2e6  # Young's modulus N/m^2
E2 = E1
nu12 = 0.0  # Poisson's ratio
nu21 = 0.0  # Poisson's ratio
G12 = E1/(2*(1+nu21)) # Shear modulus
rho = 8.0e3  # *9.807 # kg/m^3
t = 100. * mm


tagPlaneStress1 = 1
ops.nDMaterial("ElasticOrthotropicPlaneStress", tagPlaneStress1, E1, E2, nu12, nu21, G12, rho)



deg2rad = pi / 180

matTags = [1]
thickness = [t]
θ = [0 * deg2rad]

gFact = [0.0, 0.0, 0.0 * 9.807]


Nlayers = len(θ)

patchTag=1
nodeStartTag=1


# Flipping control point to u,v format
controlPts = surf.ctrlpts2d[:]
controlPts = np.array(compatibility.flip_ctrlpts2d(controlPts))

for dim in controlPts: # Unweighting control pts
    for point in dim:
        point[0:3]/=point[3]

# Creating a Patch in OpenSees
ops.IGA("Patch", patchTag, nodeStartTag, surf.degree_u, surf.degree_v, surf.ctrlpts_size_u, surf.ctrlpts_size_v,
        "-type", "KLShell",
        # "-nonLinearGeometry", 0,
        "-planeStressMatTags", *matTags,
        "-gFact", *gFact,
        "-theta", *θ,
        "-thickness", *thickness,
        "-uKnot", *surf.knotvector_u, "-vKnot", *surf.knotvector_v, "-controlPts", *controlPts.flatten())




# Creating constraints
for n in ops.getNodeTags():
    n = int(n)
    if n in [1, 2, 13, 14]:
        ops.fix(n, 1, 1, 1)


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

# Pz = 4.0*1e-12
Pz = 4.0
for n in [12, 24]:
    # # Vertical load
    ops.load(n, 0, 0, Pz / 2.0)

I = Lb * (sum(thickness)**3) / 12.0

print("Finished loading nodes")


print("Starting analysis")

# Create test
ops.test("NormDispIncr", 1.0e-8, 500, 1)


# create SOE
ops.system("UmfPack")

# create DOF number
ops.numberer("Plain")

# create constraint handler
ops.constraints("Plain")
# ops.constraints("Penalty",1,1)


# ops.algorithm("Linear")
# ops.algorithm("Newton")
# ops.algorithm("SecantNewton")
ops.algorithm("NewtonLineSearch")
# ops.algorithm("ModifiedNewton")
# ops.algorithm("KrylovNewton")
# ops.algorithm("BFGS")
# ops.algorithm("Broyden")

# create integrator

# nSteps = 50
# ops.integrator("LoadControl", 1.0 / nSteps)

delta = 0.67/2
defMax = 6.7
nSteps = abs(int(defMax / delta))
# ops.integrator("LoadControl", 1.0 / nSteps)
ops.integrator("DisplacementControl", 12, 3, delta)


# create analysis object
ops.analysis("Static")

import matplotlib.pyplot as plt
data = np.zeros((nSteps + 1, 3))

# perform the analysis
for j in range(nSteps):
    print("=================================")
    print(f"Load step {j}")
    print("=================================")
    result = ops.analyze(1)
    if result != 0:
        break
        exit(-1)
    data[j + 1, 0] = abs(ops.nodeDisp(12, 3))
    data[j + 1, 1] = abs(ops.nodeDisp(12, 1))
    data[j + 1, 2] = abs(ops.getLoadFactor(1) * Pz)
    print(data[j + 1, 2])





# Adding deformation to control points
fdef = 1e0
nodes=ops.getNodeTags()
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



# Visualize surface
surfVisualize(surf, hold=False)


result = np.array(ops.nodeDisp(24, 3))
print("result: ", result)

F=Pz*np.arange(0,1.05,0.05)
U=[0,0.026,0.103,0.224,0.381,0.563,0.763,0.971,1.184,1.396,1.604,1.807,2.002,2.19,2.37,2.541,2.705,2.861,3.01,3.151,3.286]
W=[0,0.663,1.309,1.922,2.493,3.015,3.488,3.912,4.292,4.631,4.933,5.202,5.444,5.66,5.855,6.031,6.19,6.335,6.467,6.588,6.698]
plt.plot(data[:, 0], data[:, 2], 'or', mfc="None", label="W disp (IGA)") #w
plt.plot(W, F, '-r',label="W disp (reference)") #w

plt.plot(data[:, 1], data[:, 2], 'ob', mfc="None", label="U disp (GNL)") #u
plt.plot(U, F, '-b',label="U disp (reference)") #u

plt.xlabel('Horizontal Displacement')
plt.ylabel('Vertical Load')
plt.legend()
plt.show()

print("Done")
