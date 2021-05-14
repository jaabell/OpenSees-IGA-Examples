#  IGA CANTILEVER PLATE. Fundamental period of 1s


import numpy as np
import opensees as ops
from math import *

# Geomgl utilities for visualization and surface manipulation
from geomdl import NURBS, compatibility, operations
from surfVisualize import *


def getCtrlPtsAndWeights(surf):
    # Flip control points to [u][v] (they are given by the getters in [v][u])
    noPtsX = surf.ctrlpts_size_u
    noPtsY = surf.ctrlpts_size_v
    weights = surf.weights
    controlPts = surf.ctrlpts2d[:]
    controlPts = compatibility.flip_ctrlpts2d(controlPts)
    controlPts = (np.array(controlPts).reshape(noPtsX * noPtsY, 4)).tolist()
    noCtrPts = len(controlPts)

    # Get separate array of control points in u,v and weights

    for i in range(len(controlPts)):
        pt = controlPts[i][:]
        wt = pt[-1]
        pt[0:3] = np.array(pt[0:3]) / wt
        controlPts[i] = pt[0:3]
        weights[i] = wt

    return controlPts, weights


def generateKnotVector(deg, nPts):

    import numpy as np
    knotVector = np.zeros(nPts + deg + 1)
    nMiddle = len(knotVector) - 2 * (deg + 1)
    step = 1.0 / (nMiddle + 1)

    ini = np.zeros(deg + 1)
    if step == 0.5:
        middle = np.array([step])
    else:
        middle = np.arange(0 + step, 1 - step, step)
        middle = np.linspace(0 + step, 1 - step, nPts +
                             deg + 1 - 2 * (deg + 1))
    fin = ini + 1

    knotVector = np.copy(ini)
    knotVector = np.append(knotVector, middle)
    knotVector = np.append(knotVector, fin)

    return knotVector


La = 8  	#
Lb = 1	    #
mm = 1.0 / 1000.  # m

ops.wipe()
ops.model('basic', '-ndm', 3, '-ndf', 3)

# These are given in v,u
controlPts = np.array([
  [0.00000 *La/10 , 0.0 *Lb , 0.0 , 1.0] ,
  [0.00000 *La/10 , 1.0 *Lb , 0.0 , 1.0] ,
  [0.31250 *La/10 , 0.0 *Lb , 0.0 , 1.0] ,
  [0.31250 *La/10 , 1.0 *Lb , 0.0 , 1.0] ,
  [0.93750 *La/10 , 0.0 *Lb , 0.0 , 1.0] ,
  [0.93750 *La/10 , 1.0 *Lb , 0.0 , 1.0] ,
  [1.87500 *La/10 , 0.0 *Lb , 0.0 , 1.0] ,
  [1.87500 *La/10 , 1.0 *Lb , 0.0 , 1.0] ,
  [3.12500 *La/10 , 0.0 *Lb , 0.0 , 1.0] ,
  [3.12500 *La/10 , 1.0 *Lb , 0.0 , 1.0] ,
  [4.37500 *La/10 , 0.0 *Lb , 0.0 , 1.0] ,
  [4.37500 *La/10 , 1.0 *Lb , 0.0 , 1.0] ,
  [5.62500 *La/10 , 0.0 *Lb , 0.0 , 1.0] ,
  [5.62500 *La/10 , 1.0 *Lb , 0.0 , 1.0] ,
  [6.87500 *La/10 , 0.0 *Lb , 0.0 , 1.0] ,
  [6.87500 *La/10 , 1.0 *Lb , 0.0 , 1.0] ,
  [8.12500 *La/10 , 0.0 *Lb , 0.0 , 1.0] ,
  [8.12500 *La/10 , 1.0 *Lb , 0.0 , 1.0] ,
  [9.06250 *La/10 , 0.0 *Lb , 0.0 , 1.0] ,
  [9.06250 *La/10 , 1.0 *Lb , 0.0 , 1.0] ,
  [9.68750 *La/10 , 0.0 *Lb , 0.0 , 1.0] ,
  [9.68750 *La/10 , 1.0 *Lb , 0.0 , 1.0] ,
  [10.00000*La/10 , 0.0 *Lb , 0.0 , 1.0] ,
  [10.00000*La/10 , 1.0 *Lb, 0.0 , 1.0]
])



patchTag = 1
nodeStartTag =1
P = 5
Q = 1

# Create a BSpline surface instance
surf = NURBS.Surface()

# Set surface degrees
surf.degree_u = P
surf.degree_v = Q

# Setting control points for surface
surf.set_ctrlpts(controlPts.tolist(), 12, 2)

# Set knot vectors
uKnot = generateKnotVector(surf.degree_u, surf.ctrlpts_size_u)
vKnot = generateKnotVector(surf.degree_v, surf.ctrlpts_size_v)

surf.knotvector_u = uKnot
surf.knotvector_v = vKnot


noPtsX = surf.ctrlpts_size_u
noPtsY = surf.ctrlpts_size_v

# Visualize surface
surfVisualize(surf, hold=True)
# exit()


# nDMaterial ElasticIsotropic $nDtag_elastic $elasticidad_probeta
# $poisson_probeta
GPa=1e9
E1 = 5*GPa  # Young's modulus N/m^2
E2 = E1
nu = 0.0  # Poisson's ratio
rho = 8e3  # *9.807 # kg/m^3
t = 0.5



tagNDmat1 = 1
ops.nDMaterial("ElasticIsotropic", tagNDmat1, E1, nu, rho)

tagNDmat2 = 2
ops.nDMaterial("ElasticIsotropic", tagNDmat2, E2, nu, rho)


# nDMaterial PlateFiber $nDtag_platefiber $nDtag_elastic
tagPlaneStress1 = 3
ops.nDMaterial("PlaneStress", tagPlaneStress1, tagNDmat1)

tagPlaneStress2 = 4
ops.nDMaterial("PlaneStress", tagPlaneStress2, tagNDmat2)

deg2rad = pi / 180


matTags = [3]
thickness = [t]
θ = [0 * deg2rad]
g=9.807

I=(Lb*(sum(thickness)**3))/12.0
W=rho*g*(sum(thickness)*Lb)
elasticSolution = abs(W*La**4/(8*E1*I))


gFact = [0.0, 0.0, 0*g]


Nlayers = len(θ)

controlPts = surf.ctrlpts2d[:]
controlPts = np.array(compatibility.flip_ctrlpts2d(controlPts))


ops.IGA("Patch", patchTag, nodeStartTag, P, Q, noPtsX, noPtsY,
        "-type", "KLShell",
        # "-nonLinearGeometry", 0,
        "-planeStressMatTags", *matTags,
        "-gFact", *gFact,
        "-theta", *θ,
        "-thickness", *thickness,
        "-uKnot", *uKnot, "-vKnot", *vKnot, "-controlPts", *controlPts.flatten())

# exit()


for n in ops.getNodeTags():
    if n in [1,2,13,14]:
        ops.fix(n,1,1,1)
    else:
        ops.fix(n,0,1,0)

# # #Fijar nodos 1, 2, 3, 4
# for n in [1,2,3,4]:
#     ops.fix(n,1,1,1)

# print("\n\n\nPRINTING DOMAIN-----------------------")
# ops.printModel()
# print("\n\n\nDONE PRINTING DOMAIN-----------------------")

# ------------------------------
# Start of analysis generation
# ------------------------------

# create TimeSeries
# ops.timeSeries("Constant", 1)
ops.timeSeries("Linear", 1)





# create a plain load pattern
ops.pattern("Plain", 1, 1 )

# Pz=3000
# for n in [12,24]:
#     ops.load(n,0,0,Pz)

weight = [0.0, 0.0, -g]
ops.eleLoad("-ele", 1, "-type", "-SelfWeight", *weight)


ops.test("EnergyIncr", 1.0e-12, 30, 0) 
ops.algorithm("Newton")
# ops.algorithm("NewtonLineSearch")
# ops.algorithm("Linear")
ops.numberer("RCM")
ops.constraints("Plain")
ops.integrator("LoadControl", 1.0,1,1.0,10.0)
ops.system("UmfPack")
ops.analysis("Static")



# perform the analysis
ops.analyze(10)

print(ops.nodeDisp(12,3))
# exit()

   
# Plot deformed surface
controlPts = surf.ctrlpts2d[:]
controlPts = compatibility.flip_ctrlpts2d(controlPts)

fdef = 1e0
i = 1
for dim in controlPts:
    for point in dim:
        point[0] += fdef * ops.nodeDisp(i,1)
        point[1] += fdef * ops.nodeDisp(i,2)
        point[2] += fdef * ops.nodeDisp(i,3)
        i+=1

nPoints = surf.ctrlpts_size_u * surf.ctrlpts_size_v
shape = np.array(compatibility.flip_ctrlpts2d(surf.ctrlpts2d[:])).shape
controlPts = np.array(controlPts).reshape(shape)
controlPts = np.array(compatibility.flip_ctrlpts2d(controlPts))

surf.set_ctrlpts(controlPts.reshape(nPoints, 4).tolist(), surf.ctrlpts_size_u, surf.ctrlpts_size_v)

# Visualize surface
surfVisualize(surf, hold=False)
# exit()





# New analysis

ops.wipeAnalysis()
ops.setTime(0.0)
ops.remove("loadPattern",1)

ops.test("EnergyIncr", 1.0e-10, 80, 0) 
# ops.test("NormUnbalance", 1.0e-10, 90, 0) 
# ops.test("NormDispIncr", 1.0e-10, 20, 0) 
ops.algorithm("Newton")
# ops.algorithm("NewtonLineSearch")
# ops.algorithm("NewtonLineSearch", 'type', 'Bisection')
ops.numberer("RCM")
ops.constraints("Plain")
ops.integrator("Newmark", 0.5, 0.25)
ops.system("UmfPack")
ops.analysis("Transient")

# Create recorder
ops.recorder('Node', '-file', 'Node2.out', '-closeOnWrite','-time', '-node', *[12], '-dof', *[3], *['disp'])
# ops.recorder('Node', '-file', 'Node.out', '-time', '-node', *[12], '-dof', *[3], *['disp'])



tFinal=5 #seconds
deltaT=0.025
nSteps=int(tFinal/deltaT)
ops.analyze(nSteps ,deltaT)



controlPts = surf.ctrlpts2d[:]
controlPts = compatibility.flip_ctrlpts2d(controlPts)

i = 1
for dim in controlPts:
    for point in dim:
        point[0] = ops.nodeCoord(i,1) + fdef * ops.nodeDisp(i,1)
        point[1] = ops.nodeCoord(i,2) + fdef * ops.nodeDisp(i,2)
        point[2] = ops.nodeCoord(i,3) + fdef * ops.nodeDisp(i,3)
        i+=1

nPoints = surf.ctrlpts_size_u * surf.ctrlpts_size_v
shape = np.array(compatibility.flip_ctrlpts2d(surf.ctrlpts2d[:])).shape
controlPts = np.array(controlPts).reshape(shape)
controlPts = np.array(compatibility.flip_ctrlpts2d(controlPts))

# print(controlPts)

surf.set_ctrlpts(controlPts.reshape(nPoints, 4).tolist(), surf.ctrlpts_size_u, surf.ctrlpts_size_v)

# Visualize surface
surfVisualize(surf, hold=False)

print(ops.nodeDisp(8,3))

# ops.record()
print("Done")

