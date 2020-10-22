import numpy as np
import opensees as ops
from math import *

# Geomgl utilities for visualization and surface manipulation
from geomdl import NURBS, compatibility
from surfVisualize import *


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


ops.model('basic', '-ndm', 3, '-ndf', 3)


uKnot = np.array([0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, ])
vKnot = np.array([0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, ])
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


# Visualize surface
# surfVisualize(surf, hold=True)

noPtsX = surf.ctrlpts_size_u
noPtsY = surf.ctrlpts_size_v


# nDMaterial ElasticIsotropic $nDtag_elastic $elasticidad_probeta
# $poisson_probeta
E = 2.03e11  # Young's modulus N/m^2
nu = 0.3  # Poisson's ratio
rho = 7.7e03  # *9.807 # kg/m^3
t = 0.05


tagNDmat1 = 1
ops.nDMaterial("ElasticIsotropic", tagNDmat, E1, nu, rho)

tagNDmat2 = 2
ops.nDMaterial("ElasticIsotropic", tagNDmat, E2, nu, rho)


# nDMaterial PlateFiber $nDtag_platefiber $nDtag_elastic
tagPlaneStress1 = 3
ops.nDMaterial("PlaneStress", tagPlaneStress, tagNDmat)

tagPlaneStress2 = 4
ops.nDMaterial("PlaneStress", tagPlaneStress, tagNDmat)

matTags =   [   3,    4,     3,     4,     3,     4,     3]
thickness = [1*mm, 2*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm]
θ =         [  0 ,   45,    90 ,  -45,     0,    45,   90]


Nlayers = len(thickness) = len(θ)

controlPts = surf.ctrlpts2d[:]
controlPts = np.array(compatibility.flip_ctrlpts2d(controlPts))


ops.IGA("Patch", patchTag, P, Q, noPtsX, noPtsY, 
    "-type", "KLShell", 
    "-planeStressMatTags", matTags, 
    "-theta", θ, 
    "-thickness", thickness, 
    "-uKnot", *uKnot, "-vKnot", *vKnot, "-controlPts", *controlPts.flatten())

# #Fijar nodos 1, 2, 3, 4
# for n in [1,2,3,4]:
#     ops.fix(n,1,1,1)


print("\n\n\nPRINTING DOMAIN-----------------------")
ops.printModel()
print("\n\n\nDONE PRINTING DOMAIN-----------------------")

# ops.system("BandSPD")
ops.system("FullGeneral")

ops.numberer("RCM")
ops.constraints("Plain")
# ops.integrator("Newmark", 0.5, 0.25)
ops.algorithm("Linear")
ops.analysis("Transient")

# #Stiffness
ops.integrator('GimmeMCK',0.0,0.0,1.0)
ops.analyze(1,0.0)
K=ops.printA('-ret')
K=np.array(K)
N=ops.systemSize()
K.shape=(N,N)
print("K: ", K)

#Mass
ops.integrator('GimmeMCK',1.0,0.0,0.0)
ops.analyze(1,0.0)
M=ops.printA('-ret')
M=np.array(M)
M.shape=(N,N)
print("M: ", M)

from scipy.sparse.linalg import eigsh
import scipy as sp
import numpy as np
W, phi_full = eigsh(K, M=M, k=24, sigma=0, maxiter=1000000000)
W = sp.sqrt(W)
order = np.argsort(W)
W = W[order].real
T = 2 * sp.pi / W
W = 1 / T
for i in range(10):
    print("Mode: ", i)
    print("W[i]: ", W[i])
    print("T[i]: ", T[i])
    print("\n")

ops.system("BandSPD")
ops.integrator("Newmark", 0.5, 0.25)

exit()

# np.save("K_ops.npy",K)

# timeSeries Path $tsTag -time $times -values $vals
nodes = ops.getNodeTags()

Nnodes = len(nodes)
Neigenvalues = 24  # arpack can only compute N-1 eigvals

w2s = ops.eigen(Neigenvalues)
order = np.argsort(w2s)
w2s = np.array(w2s, dtype=np.float64)[order]
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


# ops.printModel()

# pattern Plain 1 $tsTag {
#     # source "Prob4_base.pullnodes.tcl"
#     # source "Prob4_base.pullnodes_load.tcl"
#     load $mothernode 0 0 0.25 0 0 0   ;# El 0.25 aqui es porque es 1/4 de probeta
# }
# constraints Transformation
# numberer RCM
# system UmfPack
# test NormDispIncr 1.0e-6 20 2
# algorithm Newton
# # algorithm NewtonLineSearch -type Secant
# # algorithm NewtonLineSearch -type Bisection
# # algorithm KrylovNewton
# algorithm BFGS
# integrator LoadControl [expr 1./$Nsteps]
# analysis Static
