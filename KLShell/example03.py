import numpy as np
import opensees as ops
from math import *

# Geomgl utilities for visualization and surface manipulation
from geomdl import NURBS, compatibility
from geomdl import operations
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
Lb = 10.0  	#
t = 0.05  	# m
mm = 1.0 / 1000.  # m


# remove existing model
ops.wipe()

# set modelbuilder
ops.model('basic', '-ndm', 3, '-ndf', 3)


uKnot = np.array([0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0])
vKnot = np.array([0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0])
controlPts = np.array([
    [-La / 2, -Lb / 2, 0, 1],
    [-La / 2, -Lb / 2 / 2, 0, 1],
    [-La / 2, Lb / 2 / 2, 0, 1],
    [-La / 2, Lb / 2, 0, 1],
    [-La / 2 / 2, -Lb / 2, 0, 1],
    [-La / 2 / 2, -Lb / 2 / 2, 0, 1],
    [-La / 2 / 2, Lb / 2 / 2, 0, 1],
    [-La / 2 / 2, Lb / 2, 0, 1],
    [La / 2 / 2, -Lb / 2, 0, 1],
    [La / 2 / 2, -Lb / 2 / 2, 0, 1],
    [La / 2 / 2, Lb / 2 / 2, 0, 1],
    [La / 2 / 2, Lb / 2, 0, 1],
    [La / 2, -Lb / 2, 0, 1],
    [La / 2, -Lb / 2 / 2, 0, 1],
    [La / 2, Lb / 2 / 2, 0, 1],
    [La / 2, Lb / 2, 0, 1]
])

# print("controlPts.flatten(): ", controlPts.flatten())

patchTag = 1
P = 2
Q = 2


# Create a BSpline surface instance
surf = NURBS.Surface()

# Set surface degrees
surf.degree_u = P
surf.degree_v = Q

# Setting control points for surface
surf.set_ctrlpts(controlPts.tolist(), 4, 4)

# Set knot vectors
uKnot = generateKnotVector(surf.degree_u, 4)
vKnot = generateKnotVector(surf.degree_u, 4)
surf.knotvector_u = uKnot
surf.knotvector_v = vKnot

# # Visualize surface
# surfVisualize(surf, hold=True)

# Refine knot vectors and update geometry (refine both directions)
# u,v
refU = 0
refV = 0
operations.refine_knotvector(surf, [refU, refV])

# Get new control points and weights after refining
weights = surf.weights
controlPts = surf.ctrlpts2d[:]

# Update geometry after refining
noPtsX = surf.ctrlpts_size_u
noPtsY = surf.ctrlpts_size_v

# Flip control points to [u][v] (they are given by the getters in [v][u])
controlPts, weights = getCtrlPtsAndWeights(surf)

# Set degrees
if refU == 0:
    P = 3
else:
    P = 4
if refV == 0:
    Q = 3
else:
    Q = 4


surf.degree_u = P
surf.degree_v = Q

uKnot = generateKnotVector(surf.degree_u, surf.ctrlpts_size_u)
surf.knotvector_u = uKnot
vKnot = generateKnotVector(surf.degree_v, surf.ctrlpts_size_v)
surf.knotvector_v = vKnot

noPtsX = surf.ctrlpts_size_u
noPtsY = surf.ctrlpts_size_v

noGPs = P + 1

noCtrPts = len(controlPts)
noDofs = noCtrPts * 3  # three displacement dofs per node

# # Visualize surface
# surfVisualize(surf, hold=True)

# exit()


# nDMaterial ElasticIsotropic $nDtag_elastic $elasticidad_probeta
# $poisson_probeta
E1 = 2.1e11  # Young's modulus N/m^2
E2 = E1
nu = 0.3  # Poisson's ratio
rho = 8.0e+03  # *9.807 # kg/m^3
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

gFact = [0.0, 0.0, 0*9.807]

# matTags = [3]
# thickness = [50.0 * mm]


# θ = [0 * deg2rad]


Nlayers = len(θ)

controlPts = surf.ctrlpts2d[:]
controlPts = np.array(compatibility.flip_ctrlpts2d(controlPts))


ops.IGA("Patch", patchTag, P, Q, noPtsX, noPtsY,
        "-type", "KLShell",
        # "-nonLinearGeometry", 0,
        "-planeStressMatTags", *matTags,
        "-gFact", *gFact,
        "-theta", *θ,
        "-thickness", *thickness,
        "-uKnot", *uKnot, "-vKnot", *vKnot, "-controlPts", *controlPts.flatten())

# Fijar nodos 
# for n in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 19, 28, 37, 46, 55, 64, 73]:
for n in [1,2,3,4,5,9,13]:
    ops.fix(n, 1, 1, 1)


print("\n\n\nPRINTING DOMAIN-----------------------")
ops.printModel()
print("\n\n\nDONE PRINTING DOMAIN-----------------------")

# create TimeSeries
ops.timeSeries("Linear", 1)

# create a plain load pattern
ops.pattern("Plain", 1, 1)

# Create the nodal load - command: load nodeID xForce yForce
Pz=6000.0
nodeLoad=16
ops.load(nodeLoad, 0.0, 0.0, Pz)


# ------------------------------
# Start of analysis generation
# ------------------------------

# create SOE
ops.system("FullGeneral")

# create DOF number
ops.numberer("RCM")

# create constraint handler
ops.constraints("Plain")

# create algorithm
# ops.algorithm("Linear")
ops.algorithm("Newton")

# create integrator
nSteps=20
ops.integrator("LoadControl", 1.0/nSteps)

# Create test
ops.test("NormDispIncr",1.0e-9,20)
# ops.test("NormUnbalance",1e-8,10)

# create analysis object
ops.analysis("Static")

# perform the analysis
import matplotlib.pyplot as plt
data=np.zeros((nSteps+1,2))
for j in range(nSteps):
    ops.analyze(1)
    data[j+1,0] = ops.nodeDisp(nodeLoad,3)
    data[j+1,1] = ops.getLoadFactor(1)*Pz

plt.plot(data[:,0], data[:,1])
plt.xlabel('Vertical Displacement')
plt.ylabel('Vertical Load')
plt.show()


fDef = 1
i = 1
for dim in controlPts:
    for point in dim:
        point[:3] += fDef * np.array(ops.nodeDisp(i))
        i += 1





# controlPts = surf.ctrlpts2d[:]
# controlPts = compatibility.flip_ctrlpts2d(surf.ctrlpts2d[:])
controlPts = (np.array(controlPts).reshape(
    surf.ctrlpts_size_u * surf.ctrlpts_size_v, 4)).tolist()

# Setting control points for surface
surf.set_ctrlpts(controlPts, surf.ctrlpts_size_v, surf.ctrlpts_size_u)

# Visualize surface
surfVisualize(surf, hold=True)
print("Finished")

exit()

# # ops.system("BandSPD")
# ops.system("FullGeneral")

# ops.numberer("RCM")
# ops.constraints("Plain")
# # ops.integrator("Newmark", 0.5, 0.25)
# ops.algorithm("Linear")
# ops.analysis("Transient")

# # Stiffness
# ops.integrator('GimmeMCK', 0.0, 0.0, 1.0)
# ops.analyze(1, 0.0)
# K = ops.printA('-ret')
# K = np.array(K)
# N = ops.systemSize()
# K.shape = (N, N)
# print("K: ", K)

# # Mass
# ops.integrator('GimmeMCK', 1.0, 0.0, 0.0)
# ops.analyze(1, 0.0)
# M = ops.printA('-ret')
# M = np.array(M)
# M.shape = (N, N)
# print("M: ", M)

# from scipy.sparse.linalg import eigsh
# import scipy as sp
# import numpy as np
# W, phi_full = eigsh(K, M=M, k=24, sigma=0, maxiter=1000000000)
# W = sp.sqrt(W)
# order = np.argsort(W)
# W = W[order].real
# T = 2 * sp.pi / W
# W = 1 / T
# for i in range(10):
#     print("Mode: ", i)
#     print("W[i]: ", W[i])
#     print("T[i]: ", T[i])
#     print("\n")

# ops.system("BandSPD")
# ops.integrator("Newmark", 0.5, 0.25)

# # exit()

# # np.save("K_ops.npy",K)

# # timeSeries Path $tsTag -time $times -values $vals
# nodes = ops.getNodeTags()

# Nnodes = len(nodes)
# Neigenvalues = 30  # arpack can only compute N-1 eigvals

# w2s = ops.eigen(Neigenvalues)
# order = np.argsort(w2s)
# w2s = np.array(w2s, dtype=np.float64)[order]
# # exit()


# for i, w2 in enumerate(w2s):
#     w = sqrt(abs(w2))
#     f = w / 2 / pi
#     T = 1 / f
#     print(f"{i} {w2} {w} {f} {T} ")


# phi = np.zeros((3 * Nnodes, Neigenvalues))

# for i, n in enumerate(nodes):
#     for j in range(Neigenvalues):
#         phi[3 * i:3 * i + 3, j] = ops.nodeEigenvector(n, j + 1)


# print(f"ϕ = {phi}")


# # ops.printModel()

# # pattern Plain 1 $tsTag {
# #     # source "Prob4_base.pullnodes.tcl"
# #     # source "Prob4_base.pullnodes_load.tcl"
# #     load $mothernode 0 0 0.25 0 0 0   ;# El 0.25 aqui es porque es 1/4 de probeta
# # }
# # constraints Transformation
# # numberer RCM
# # system UmfPack
# # test NormDispIncr 1.0e-6 20 2
# # algorithm Newton
# # # algorithm NewtonLineSearch -type Secant
# # # algorithm NewtonLineSearch -type Bisection
# # # algorithm KrylovNewton
# # algorithm BFGS
# # integrator LoadControl [expr 1./$Nsteps]
# # analysis Static
