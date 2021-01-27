#  IGA Scordelis Lo Roof
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



La = 10.0  	#
Lb = 1.0  	#
t = 0.05  	# m
mm = 1.0 / 1000.  # m


ops.model('basic', '-ndm', 3, '-ndf', 3)

R=25
L=50
phi=40
deg2rad = pi / 180

controlPts=np.zeros((4,3,3))

# First zKnot

x1=np.sqrt(R**2+R**2*np.tan(deg2rad*20)**2)

point0=[R*np.cos(deg2rad*(90-phi)), R*np.sin(deg2rad*(90-phi)), 0, 1]
point1=[x1*np.cos(deg2rad*(90-phi/2)), x1*np.sin(deg2rad*(90-phi/2)), 0, 1]
point2=[0,R,0,1]
point3=[x1*np.cos(deg2rad*(90+phi/2)), x1*np.sin(deg2rad*(90+phi/2)), 0, 1]
point4=[R*np.cos(deg2rad*(90+phi)), R*np.sin(deg2rad*(90+phi)), 0, 1]


# Second zKnot
z=L/2

# point3=[R*np.cos(deg2rad*(90-phi)), R*np.sin(deg2rad*(90-phi)),z,1]
# point4=[x1*np.cos(deg2rad*70), x1*np.sin(deg2rad*70), z, 1]
# point5=[0,R,z,1]

point5=[R*np.cos(deg2rad*(90-phi)), R*np.sin(deg2rad*(90-phi)), z, 1]
point6=[x1*np.cos(deg2rad*(90-phi/2)), x1*np.sin(deg2rad*(90-phi/2)), z, 1]
point7=[0,R,z,1]
point8=[x1*np.cos(deg2rad*(90+phi/2)), x1*np.sin(deg2rad*(90+phi/2)), z, 1]
point9=[R*np.cos(deg2rad*(90+phi)), R*np.sin(deg2rad*(90+phi)), z, 1]

# Third zKnot
z=L

# point6=[R*np.cos(deg2rad*(90-phi)), R*np.sin(deg2rad*(90-phi)),z,1]
# point7=[x1*np.cos(deg2rad*70), x1*np.sin(deg2rad*70), z, 1]
# point8=[0,R,z,1]

point10=[R*np.cos(deg2rad*(90-phi)), R*np.sin(deg2rad*(90-phi)), z, 1]
point11=[x1*np.cos(deg2rad*(90-phi/2)), x1*np.sin(deg2rad*(90-phi/2)), z, 1]
point12=[0,R,z,1]
point13=[x1*np.cos(deg2rad*(90+phi/2)), x1*np.sin(deg2rad*(90+phi/2)), z, 1]
point14=[R*np.cos(deg2rad*(90+phi)), R*np.sin(deg2rad*(90+phi)), z, 1]

# Weights
fac=np.cos(deg2rad*(phi/2))
point1[3]=fac
point2[3]=fac
point3[3]=fac

point6[3]=fac
point7[3]=fac
point8[3]=fac

point11[3]=fac
point12[3]=fac
point13[3]=fac

controlPts=[point0,point1,point2,point3,point4,point5,point6,point7,point8,point9,point10,point11,point12,point13,point14]

for point in controlPts:
    for i in range(3):
        point[i]*=point[3]

# Fix points 2, 7, 12
for n in [2,7,12]:
    controlPts[n][1]/=controlPts[n][3]

# These are given in v,u


patchTag = 1
P = 2
Q = 2

# Create a BSpline surface instance
surf = NURBS.Surface()

# Set surface degrees
surf.degree_u = P
surf.degree_v = Q

# Setting control points for surface
surf.set_ctrlpts(controlPts, 3, 5)

# Set knot vectors
uKnot = generateKnotVector(surf.degree_u, surf.ctrlpts_size_u)
vKnot = generateKnotVector(surf.degree_v, surf.ctrlpts_size_v)

surf.knotvector_u = uKnot
surf.knotvector_v = vKnot


noPtsX = surf.ctrlpts_size_u
noPtsY = surf.ctrlpts_size_v

# Visualize surface
surfVisualize(surf, hold=True)
# 


# nDMaterial ElasticIsotropic $nDtag_elastic $elasticidad_probeta
# $poisson_probeta
E1 = 4.32e8  # Young's modulus N/m^2
E2 = E1
nu = 0.0  # Poisson's ratio
rho = 8.0e3  # *9.807 # kg/m^3
t = 0.05


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

matTags = [3, 4, 3, 4, 3]
thickness = [10. * mm, 10. * mm, 10. * mm, 10. * mm, 10. * mm]
θ = [0 * deg2rad, 45 * deg2rad, 90 * deg2rad, -45 * deg2rad, 0 * deg2rad]

gFact = [0.0, 0.0, 0.0 * 9.807]

matTags = [3]
thickness = [50.0 * mm]
θ = [0 * deg2rad]


Nlayers = len(θ)

controlPts = surf.ctrlpts2d[:]
controlPts = np.array(compatibility.flip_ctrlpts2d(controlPts))

print("controlPts.tolist(): ", controlPts.tolist())


ops.IGA("Patch", patchTag, P, Q, noPtsX, noPtsY,
        "-type", "KLShell",
        # "-nonLinearGeometry", 0,
        "-planeStressMatTags", *matTags,
        "-gFact", *gFact,
        "-theta", *θ,
        "-thickness", *thickness,
        "-uKnot", *uKnot, "-vKnot", *vKnot, "-controlPts", *controlPts.flatten())

exit()


for n in [1,2,3,4,5,6,7,8]:
    if n in [1,2,3,4]:
        ops.fix(n,1,1,1)
    else:
        ops.fix(n,1,1,0)

# # #Fijar nodos 1, 2, 3, 4
# for n in [1,2,3,4]:
#     ops.fix(n,1,1,1)

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
# Cargar nodos 7,8
Pz=-0.5e5
for n in [7,8]:
    ops.load(n,0,0,Pz)
print("Finished loading nodes")

# create SOE
ops.system("FullGeneral")
# ops.system("BandSPD")

# create DOF number
ops.numberer("RCM")



print("Starting analysis")

# create constraint handler
ops.constraints("Plain")
# ops.constraints("Penalty", 1e12*E1,1e12*E1)

# Create test
ops.test("NormDispIncr", 1.0e-9, 50, 1)
# ops.test("NormUnbalance",1e-8,10)

# create algorithm

# ops.algorithm("Linear")
ops.algorithm("Newton")
# ops.algorithm("ModifiedNewton")
# ops.algorithm("KrylovNewton")

# create integrator
# nSteps=10
# ops.integrator("LoadControl", 1.0/nSteps)
ops.integrator("LoadControl", 1.0)


# create analysis object
ops.analysis("Static")

ops.analyze(1)

print("Finished analysis")

controlPts = surf.ctrlpts2d[:]
controlPts = compatibility.flip_ctrlpts2d(controlPts)

fDef = 1e0
i = 1
for dim in controlPts:
    for point in dim:
        point[:3] += fDef * np.array(ops.nodeDisp(i))
        i += 1

controlPts = compatibility.flip_ctrlpts2d(controlPts)

controlPts = (np.array(controlPts).reshape(
    surf.ctrlpts_size_u * surf.ctrlpts_size_v, 4)).tolist()

# Setting control points for surface
print("controlPts: ", controlPts)
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

print("ops.nodeDisp(7,2): ", 1000*np.array(ops.nodeDisp(7)),"mm")
print("ops.nodeDisp(8,2): ", 1000*np.array(ops.nodeDisp(8)),"mm")

print("Done")
