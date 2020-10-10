import numpy as np
import opensees as ops

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

P, Q = 2, 1
noPtsU, noPtsV = 3, 13
uKnot = [0, 0, 0, 1, 1, 1]
vKnot = [0, 0, 0.0833333333333, 0.166666666667, 0.25, 0.333333333333, 0.416666666667,
         0.5, 0.583333333333, 0.666666666667, 0.75, 0.833333333333, 0.916666666667, 1, 1]

controlPts = [
    [-0.056987131848, -0.007203602858,  -0.229, 1],
    [-0.057499999937, 1.913055e-006, -0.229, 1],
    [-0.05698672389, 0.007207412846, -0.229, 1],
    [-0.057047829316, -0.007273454447,  -0.119, 1],
    [-0.057499961249, -6.7861588e-005,  -0.119, 1],
    [-0.057064583187, 0.007138832831, -0.119000000001, 1],
    [-0.079794778661, 0.021369498005, -0, 1],
    [-0.079524941035, 0.033126988323, 0, 1],
    [-0.067301719459, 0.038293031987, -1e-012, 1],
    [-0.073751655992, 0.005336416859, 0.24, 1],
    [-0.075568637924, 0.015533595701, 0.24, 1],
    [-0.06569805411, 0.022105408648, 0.24, 1],
    [-0.067174927869, -0.003682525739,  0.48, 1],
    [-0.069584388864, 0.00543290258, 0.48, 1],
    [-0.061821473935, 0.012662984031, 0.48, 1],
    [-0.063374336449, -0.007710981686,  0.72, 1],
    [-0.066816460403, 0.00077594064, 0.72, 1],
    [-0.059129961601, 0.008343910706, 0.72, 1],
    [-0.051552924986, -0.006565078205,  0.96, 1],
    [-0.054605146857, 0.000683553987, 0.96, 1],
    [-0.047969353108, 0.007320587936, 0.96, 1],
    [-0.04168665773, -0.005678364952, 1.2, 1],
    [-0.044417107237, 0.00063222129, 1.2, 1],
    [-0.038556984706, 0.006615460774, 1.2, 1],
    [-0.028901292344, -0.004692571975,  1.44, 1],
    [-0.031236232726, 0.000559037969, 1.44, 1],
    [-0.026295856415, 0.005748733282, 1.44, 1],
    [-0.019331913844, -0.003786905972,  1.68, 1],
    [-0.021246484266, 0.000464496821, 1.68, 1],
    [-0.017239590984, 0.004830909624, 1.68, 1],
    [-0.014462494662, -0.002920852489,  1.92, 1],
    [-0.015957931846, 0.00031814018, 1.92, 1],
    [-0.012932450378, 0.003757703925, 1.92, 1],
    [-0.01004992451, -0.002197338409, 2.16, 1],
    [-0.011138889382, 0.000208443071, 2.16, 1],
    [-0.008935568477, 0.00281901551, 2.16, 1],
    [0.026173908381, 0.000732170416, 2.26, 1],
    [0.025956115407, 0.001213326712, 2.26, 1],
    [0.026396779588, 0.0017354412,  2.26, 1]
]

# print("controlPts.flatten(): ", controlPts.flatten())

patchTag = 1

# Create a BSpline surface instance
surf = NURBS.Surface()

# Set surface degrees
surf.degree_u = P
surf.degree_v = Q

# # Setting control points for surface
# surf.set_ctrlpts(controlPts.tolist(), noPtsU, noPtsV)

# Set control points
surf.set_ctrlpts(controlPts, noPtsV, noPtsU)

test=compatibility.flip_ctrlpts2d(surf.ctrlpts2d)
# test=surf.ctrlpts2d[:]
test=np.array(test).reshape(noPtsU*noPtsV,4)
test=test.tolist()

surf.set_ctrlpts(test, noPtsU, noPtsV)

surf.ctrlpts_size_u = noPtsU
surf.ctrlpts_size_v = noPtsV


# Set knot vectors
# uKnot = generateKnotVector(surf.degree_u, 4)
# vKnot = generateKnotVector(surf.degree_v, 4)
surf.knotvector_u = uKnot
surf.knotvector_v = vKnot


# Visualize surface
surfVisualize(surf, hold=True)

noPtsX = surf.ctrlpts_size_u
noPtsY = surf.ctrlpts_size_v

controlPts=np.array(controlPts)

# ops.IGA("Patch", patchTag, P, Q, "-uKnot", *uKnot, "-vKnot", *vKnot,
# "-controlPts", *controlPts.flatten())
ops.IGA("Patch", patchTag, P, Q, noPtsU, noPtsV, "-uKnot", *uKnot, "-vKnot", *vKnot, "-controlPts", *controlPts.flatten())
# print("noPtsX: ", noPtsX)