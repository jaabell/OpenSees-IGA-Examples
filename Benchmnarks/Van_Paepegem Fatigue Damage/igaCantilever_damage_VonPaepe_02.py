#  IGA CANTILEVER ROTATING


import matplotlib
from matplotlib.pylab import *
import numpy as np
import opensees as ops
from math import *

# Geomgl utilities for visualization and surface manipulation
from geomdl import NURBS, compatibility, operations, knotvector, multi
from surfVisualize import *

mm = 1.0/1000

Ltotal = 154 * mm
L = 54*mm
a = (Ltotal-L)/1.15

b = 30 * mm

# Creating first surface (from clamp to end)

# These are given in v,u
controlPts_L = np.array([
    [0 * L, 0.0 * b, 0.0, 1.0],
    [0 * L, 1.0 * b, 0.0, 1.0],
    [1.0/5.0 * L, 0.0 * b, 0.0, 1.0],
    [1.0/5.0 * L, 1.0 * b, 0.0, 1.0],
    [2.0/5.0 * L, 0.0 * b, 0.0, 1.0],
    [2.0/5.0 * L, 1.0 * b, 0.0, 1.0],
    [3.0/5.0 * L, 0.0 * b, 0.0, 1.0],
    [3.0/5.0 * L, 1.0 * b, 0.0, 1.0],
    [4.0/5.0 * L, 0.0 * b, 0.0, 1.0],
    [4.0/5.0 * L, 1.0 * b, 0.0, 1.0],
    [5.0/5.0 * L, 0.0 * b, 0.0, 1.0],
    [5.0/5.0 * L, 1.0 * b, 0.0, 1.0],
    [L+a/2, 0.0 * b, 0.0, 1.0],
    [L+a/2, 1.0 * b, 0.0, 1.0],
    [L+a, 0.0 * b, 0.0, 1.0],
    [L+a, 1.0 * b, 0.0, 1.0]

])

# Create a Nurbs surface instance
surf_L = NURBS.Surface()

# Set surface degrees
surf_L.degree_u = 2
surf_L.degree_v = 1

# Setting control points for surface
surf_L.set_ctrlpts(controlPts_L.tolist(), 8, 2)

# Set knot vectors
surf_L.knotvector_u = knotvector.generate(
    surf_L.degree_u, surf_L.ctrlpts_size_u)
surf_L.knotvector_v = knotvector.generate(
    surf_L.degree_v, surf_L.ctrlpts_size_v)


# Creating second surface (from clamp to end)

# These are given in v,u
controlPts_La = np.array([
    [L, 0, 0, 1],
    [L, b, 0, 1],
    [L+a/2, 0, 0, 1],
    [L+a/2, b, 0, 1],
    [L+a, 0, 0, 1],
    [L+a, b, 0, 1]
])

# Create a Nurbs surface instance
surf_La = NURBS.Surface()

# Set surface degrees
surf_La.degree_u = 2
surf_La.degree_v = 1

# Setting control points for surface
surf_La.set_ctrlpts(controlPts_La.tolist(), 3, 2)

# Set knot vectors
surf_La.knotvector_u = knotvector.generate(
    surf_La.degree_u, surf_La.ctrlpts_size_u)
surf_La.knotvector_v = knotvector.generate(
    surf_La.degree_v, surf_La.ctrlpts_size_v)


# Creating container for multipatches
surfList = [surf_L, surf_La]

container = multi.SurfaceContainer(surfList)

# Visualize surface

container.sample_size = 30
for surf in container:
    surf.evaluate()

# Visualization configuration
container.vis = VisVTK.VisSurface(ctrlpts=True, legend=False)
# container.vis.ctrlpts_offset=0.1

# Render the surface
evalcolor = ["green", 'blue']
container.render(evalcolor=evalcolor)


# Create OpenSees model

ops.wipe()
ops.model('basic', '-ndm', 3, '-ndf', 3)


# Material parameters
GPa = 1e9
E1 = 24.57*GPa  # Young's modulus N/m^2
E2 = 23.94*GPa
G12 = 4.83*GPa  # Shear modulus
nu12 = 0.153  # Poisson's ratio
nu21 = E2*nu12/E1  # Poisson's ratio
rho = 1e3  # kg/m^3, won't use it

MPa = 1e6
Xt = 390.7*MPa
Xc = 345.1*MPa
Yt = 390.7*MPa
Yc = 345.1*MPa
S = 100.6*MPa
c1 = 0.003
c2 = 30.0
c3 = 3.5e-6
c4 = 0.85
c5 = 93
c6 = 0
c7 = 0
c8 = 0
c9 = 0.6
b = 1.0

tagPlaneStress1 = 1
vonPaepeParams = [E1, E2, nu12, nu21, G12, rho, Xt, Xc,
                  Yt, Yc, S, c1, c2, c3, c4, c5, c6, c7, c8, c9, b]
ops.nDMaterial("VonPapaDamage", tagPlaneStress1, *vonPaepeParams)


tagNDmat2 = 2
ops.nDMaterial("ElasticIsotropic", tagNDmat2, 100*E1, nu12, rho)
tagPlaneStress2 = 3
ops.nDMaterial("PlaneStress", tagPlaneStress2, tagNDmat2)

materialTags = [tagPlaneStress1, tagPlaneStress2]


# Laminated parameters

deg2rad = pi / 180

nLayers = 8
totalThickness = 2.72*mm
thick_i = totalThickness/nLayers
thickness = [thick_i]*nLayers
θ = deg2rad*np.array([0]*nLayers)
matTags = [materialTags[0]]*nLayers


names = ["surf_L", "surf_La"]
shellType = 'KLShell'
patchTag = 1
nodeStartTag = 1

nodesMap = []

for i in range(len(container)):
    surf = container[i]
    name = names[i]

    print("Creating patch ", name)

    if name == 'surf_La':
        matTags = [materialTags[1]]
        θ = [0]
        thickness = [totalThickness]

    # Flipping control point to u,v format
    controlPts = surf.ctrlpts2d[:]
    controlPts = np.array(compatibility.flip_ctrlpts2d(controlPts))

    for dim in controlPts:  # Unweighting control pts
        for point in dim:
            point[0:3] /= point[3]


    print(patchTag)
    print(nodeStartTag)
    print(surf.degree_u)
    print(surf.degree_v)
    print(surf.ctrlpts_size_u)
    print(surf.ctrlpts_size_v)
    print(shellType)
    print(matTags)
    print(θ)
    print(thickness)
    print(surf.knotvector_u)
    print(surf.knotvector_v)

    # Creating a Patch in OpenSees
    ops.IGA("Patch", patchTag, nodeStartTag, surf.degree_u, surf.degree_v, surf.ctrlpts_size_u, surf.ctrlpts_size_v,
            "-type", shellType,
            # "-nonLinearGeometry", 0,
            "-planeStressMatTags", *matTags,
            "-theta", *θ,
            "-thickness", *thickness,
            "-uKnot", *surf.knotvector_u, "-vKnot", *surf.knotvector_v, "-controlPts", *controlPts.flatten())

    # Get the nodes on current patch
    nodesMap.append(np.arange(nodeStartTag, ops.getNodeTags()[-1] + 1).tolist())

    # Update patchTag, nodeStartTag
    lastElTag = ops.getEleTags()[-1]
    lastNodeTag = ops.getNodeTags()[-1]
    patchTag = lastElTag + 1
    nodeStartTag = lastNodeTag + 1


# Creating DOF numberer
ops.numberer("RCM")

# Creating constraints handler
ops.constraints("Plain")


# Creating equalDofs for patches
ops.equalDOF(8, 19, 1, 2, 3)
ops.equalDOF(7, 18, 1, 2, 3)
ops.equalDOF(6, 17, 1, 2, 3)

ops.equalDOF(16, 22, 1, 2, 3)
ops.equalDOF(15, 21, 1, 2, 3)
ops.equalDOF(14, 20, 1, 2, 3)


# Creating constraints
for n in ops.getNodeTags():
    if n in [1, 2, 9, 10]:
        ops.fix(n, 1, 1, 1)

print("\n\n\nPRINTING DOMAIN-----------------------")
ops.printModel()
print("\n\n\nDONE PRINTING DOMAIN-----------------------")

print("Succesfully created model")


ω = 2*pi  # rad/s
tMax = 1  # 5 seconds
deltaT = 2.5e-2
t = np.arange(0, tMax+deltaT, deltaT)

uMax = 30.4*mm/2
# uMax = 1.2*mm/2
uTip = np.sin(ω*t-pi/2)*uMax+uMax

plot(t, 1000*uTip, '-o')
show()


# create TimeSeries
# ops.timeSeries("Path", 1, '-time', *(t.tolist()), '-values', *(uTip.tolist()))
ops.timeSeries("Linear", 1)
# Crear time series trigonometrico de tiempo arbitrario

# create a plain load pattern
ops.pattern("Plain", 1, 1)

# Loading tip nodes
load = 0.5
ops.load(8, 0, 0, load)
ops.load(16, 0, 0, load)

# print("\n\n\nPRINTING DOMAIN-----------------------")
ops.printModel()
# print("\n\n\nDONE PRINTING DOMAIN-----------------------")


# ------------------------------
# Start of analysis generation
# ------------------------------


# # create a plain load pattern for self weight
# ops.timeSeries("Constant", 5)

# ops.pattern("Plain", 5, 5)

# g = 9.807
# weight = [0.0, 0.0, -g]
# ops.eleLoad("-ele", 1, "-type", "-SelfWeight", *weight)

# ops.system("UmfPack")
ops.system("FullGeneral")

# Analysis
# ops.test("EnergyIncr", 1.0e-7, 100, 0)
# ops.test("NormUnbalance", 1.0e-7, 100, 1)
ops.test("NormDispIncr", 1.0e-7, 100, 0)



# Creating algorithm
ops.algorithm("Newton")
# ops.algorithm("Linear")
# ops.algorithm("NewtonLineSearch")
# ops.algorithm("NewtonLineSearch", 'type', 'Bisection')


# Create recorder
# ops.recorder('Node', '-file', 'Node12_Z.out', '-closeOnWrite', '-time', '-node', *[12], '-dof', *[3], *['disp'])

# Compute cycle
def computeCycle():
    D0 = 0
    loadFactor_max = 0
    for j in range(1, len(t)):
        delta = uTip[j]-D0
        # print("delta = ", delta)

        # Creating integrator
        ops.integrator("DisplacementControl", 16, 3, delta)

        # Create analysis type
        ops.analysis("Static")

        if (ops.analyze(1) != 0):
            print("Analysis failed")
            exit()
        # exit()
        # print("disp12 = ", ops.nodeDisp(12, 3),'\n')
        loadFactor = ops.getLoadFactor(1)
        loadFactor_max = max(loadFactor, loadFactor_max)
        # print("loadFactor = ", loadFactor)
        D0 = uTip[j]
    return loadFactor_max


data = np.array([
    [0, 106.49681528662421],
    [1186.9436201780336, 104.39490445859873],
    [2373.887240356067, 102.05944798301488],
    [3560.830860534101, 100.19108280254778],
    [9495.548961424298, 98.32271762208069],
    [16617.210682492558, 97.38853503184714],
    [27299.703264094976, 96.4543524416136],
    [46290.80118694363, 95.75371549893845],
    [68842.72997032641, 95.28662420382166],
    [89020.7715133531, 95.05307855626327],
    [108011.86943620178, 95.05307855626327],
    [127002.96735905044, 94.5859872611465],
    [148367.95252225522, 94.5859872611465],
    [172106.824925816, 94.11889596602973],
    [205341.24629080118, 93.65180467091295],
    [231454.0059347181, 93.41825902335458],
    [263501.4836795252, 93.18471337579618],
    [316913.94658753707, 92.7176220806794],
    [351335.3115727003, 92.25053078556265],
    [376261.1275964392, 91.54989384288749],
    [394065.28189910983, 91.31634819532908],
    [420178.0415430267, 90.84925690021232],
    [468842.7299703264, 89.91507430997876],
    [499703.2640949555, 89.21443736730362],
    [530563.7982195846, 88.51380042462847],
    [562611.2759643918, 87.57961783439491],
    [594658.7537091989, 86.64543524416136],
    [626706.231454006, 84.77707006369428],
    [658753.7091988132, 82.90870488322719],
    [690801.1869436203, 81.50743099787687],
    [721661.7210682493, 78.47133757961785],
    [752522.2551928784, 74.03397027600849],
    [779821.9584569733, 68.89596602972401],
    [798813.0563798221, 65.3927813163482]
])

cycles_data = data[:, 0]
force_data = data[:, 1]
steps = np.diff(cycles_data)
midSteps = np.diff(np.linspace(1, 1186, 50))
steps = np.sort(np.concatenate(
    [[1]*10, [10]*10, [100]*10, [200]*10, [500]*10, [1000]*10, midSteps, steps, [1000]*10]))
steps=np.concatenate([[1]*10,[10]*10,[100]*10,[500]*10,[1000]*10])

# Computing loadFactors
nSteps = len(steps)
loadFactors = np.zeros(nSteps, dtype=np.float64)
cycles = np.zeros(nSteps)
for i in range(nSteps):
    print("\n\nStep = ", i, " of ", nSteps)
    loadFactor_max = computeCycle()
    print('loadFactor = ', loadFactor_max)
    loadFactors[i] = loadFactor_max
    nCycles = steps[i]
    print("Advancing damage state in ",nCycles,"cycles")
    # exit()
    ops.setParameter('-val', int(nCycles), "-ele", 1, "advanceDamageState")
    ops.setParameter('-val', 0, "-ele", 1, "resetMaxStress")
    if i == 0:
        cycles[i] = nCycles
    else:
        cycles[i] = cycles[i-1]+nCycles

matplotlib.rc('axes.formatter', useoffset=False)

print(loadFactors)

# plt.rcParams['axes.formatter.useoffset'] = False
plot(cycles, loadFactors, 'or')
plot(cycles_data, force_data, '-b')

# ticklabel_format(useOffset=False)
show()
print("Done")
