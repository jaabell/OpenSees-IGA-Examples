import numpy as np
import opensees as ops

# Geomgl utilities for visualization and surface manipulation
from geomdl import NURBS
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



# nDMaterial ElasticIsotropic $nDtag_elastic $elasticidad_probeta $poisson_probeta
E = 2.03e11  # Young's modulus N/m^2
nu = 0.3  # Poisson's ratio
rho = 7.7e03 #*9.807 # kg/m^3 
t = 0.0254


tagNDmat = 1
ops.nDMaterial("ElasticIsotropic", tagNDmat, E, nu, rho)

# nDMaterial PlateFiber $nDtag_platefiber $nDtag_elastic
tagPlateFiber = 2
ops.nDMaterial("PlateFiber",tagPlateFiber, tagNDmat)


# section PlateFiber $secTag_probeta $nDtag_platefiber $espesor_probeta
tagSection = 1
ops.section("PlateFiber", tagSection, tagPlateFiber, t)





ops.IGA("Patch", patchTag, P, Q, noPtsX, noPtsY, "-type", "KLShell", "-sectionTag", tagSection, "-uKnot", *uKnot, "-vKnot", *vKnot, "-controlPts", *controlPts.flatten())

print("\n\n\nPRINTING DOMAIN-----------------------")
ops.printModel()
print("\n\n\nDONE PRINTING DOMAIN-----------------------")



# timeSeries Path $tsTag -time $times -values $vals 
ops.eigen(10)

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