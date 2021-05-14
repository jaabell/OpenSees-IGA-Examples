#  IGA CANTILEVER ROTATING


from matplotlib.pylab import *
import numpy as np
import opensees as ops
from math import *

# Geomgl utilities for visualization and surface manipulation
from geomdl import NURBS, compatibility, operations, knotvector
from surfVisualize import *


mm = 1.0/1000
La = 160 * mm    #
# La = 54 * mm    #
Lb = 30 * mm

# These are given in v,u
controlPts = np.array([
    [0.00000 / 10 * La, 0.0 * Lb, 0.0, 1.0],
    [0.00000 / 10 * La, 1.0 * Lb, 0.0, 1.0],
    [0.31250 / 10 * La, 0.0 * Lb, 0.0, 1.0],
    [0.31250 / 10 * La, 1.0 * Lb, 0.0, 1.0],
    [0.93750 / 10 * La, 0.0 * Lb, 0.0, 1.0],
    [0.93750 / 10 * La, 1.0 * Lb, 0.0, 1.0],
    [1.87500 / 10 * La, 0.0 * Lb, 0.0, 1.0],
    [1.87500 / 10 * La, 1.0 * Lb, 0.0, 1.0],
    [3.12500 / 10 * La, 0.0 * Lb, 0.0, 1.0],
    [3.12500 / 10 * La, 1.0 * Lb, 0.0, 1.0],
    [4.37500 / 10 * La, 0.0 * Lb, 0.0, 1.0],
    [4.37500 / 10 * La, 1.0 * Lb, 0.0, 1.0],
    [5.62500 / 10 * La, 0.0 * Lb, 0.0, 1.0],
    [5.62500 / 10 * La, 1.0 * Lb, 0.0, 1.0],
    [6.87500 / 10 * La, 0.0 * Lb, 0.0, 1.0],
    [6.87500 / 10 * La, 1.0 * Lb, 0.0, 1.0],
    [8.12500 / 10 * La, 0.0 * Lb, 0.0, 1.0],
    [8.12500 / 10 * La, 1.0 * Lb, 0.0, 1.0],
    [9.06250 / 10 * La, 0.0 * Lb, 0.0, 1.0],
    [9.06250 / 10 * La, 1.0 * Lb, 0.0, 1.0],
    [9.68750 / 10 * La, 0.0 * Lb, 0.0, 1.0],
    [9.68750 / 10 * La, 1.0 * Lb, 0.0, 1.0],
    [10.00000 / 10 * La, 0.0 * Lb, 0.0, 1.0],
    [10.00000 / 10 * La, 1.0 * Lb, 0.0, 1.0]
])


patchTag = 1
nodeStartTag = 1
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
uKnot = knotvector.generate(surf.degree_u, surf.ctrlpts_size_u)
vKnot = knotvector.generate(surf.degree_v, surf.ctrlpts_size_v)

surf.knotvector_u = uKnot
surf.knotvector_v = vKnot


noPtsX = surf.ctrlpts_size_u
noPtsY = surf.ctrlpts_size_v

# Visualize surface
surfVisualize(surf, hold=False)



# Create OpenSees model

ops.wipe()
ops.model('basic', '-ndm', 3, '-ndf', 3)


# Material parameters
GPa  = 1e9
E1   = 24.57*GPa  # Young's modulus N/m^2
E2   = 23.94*GPa
G12  = 4.83*GPa  # Shear modulus
nu12 = 0.153  # Poisson's ratio
nu21 = E2*nu12/E1  # Poisson's ratio
rho  = 1e3  # kg/m^3, won't use it

MPa  = 1e6
Xt   = 390.7*MPa
Xc   = 345.1*MPa
Yt   = 390.7*MPa
Yc   = 345.1*MPa
S    = 100.6*MPa
c1   = 0.003
c2   = 30.0
c3   = 3.5e-6
c4   = 0.85
c5   = 93
c6   = 0
c7   = 0
c8   = 0
c9   = 0.6
b    = 1.0

tagPlaneStress = 1
vonPaepeParams = [E1, E2, nu12, nu21, G12, rho, Xt, Xc,
                  Yt, Yc, S, c1, c2, c3, c4, c5, c6, c7, c8, c9, b]
ops.nDMaterial("VonPapaDamage", tagPlaneStress, *vonPaepeParams)
print("Created Von-Paepe material")
deg2rad = pi / 180

nLayers = 8
totalThickness = 2.72*mm
thick_i = totalThickness/nLayers
thickness = [thick_i]*nLayers
θ = deg2rad*np.array([0]*nLayers)
matTags = [1]*nLayers


controlPts = surf.ctrlpts2d[:]
controlPts = np.array(compatibility.flip_ctrlpts2d(controlPts))


ops.IGA("Patch", patchTag, nodeStartTag, P, Q, noPtsX, noPtsY,
        "-type", "KLShell",
        # "-nonLinearGeometry", 0,
        "-planeStressMatTags", *matTags,
        "-theta", *θ,
        "-thickness", *thickness,
        "-uKnot", *uKnot, "-vKnot", *vKnot, "-controlPts", *controlPts.flatten())


# ops.setParameter('-val', 0, "-ele", 1, "resetMaxStress")

# ops.updateParameter()

print("Succesfully created model")

# Creating equalDofs for tip
ops.equalDOF(12, 11, 3)
ops.equalDOF(24, 23, 3)

# Creating constraints
for n in ops.getNodeTags():
    if n in [1, 2, 13, 14]:
        ops.fix(n, 1, 1, 1)

    # else:
    #     ops.fix(n, 0, 1, 0)


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
ops.load(12, 0, 0, load)
ops.load(24, 0, 0, load)

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


# Analysis
# ops.test("EnergyIncr", 1.0e-7, 100, 0)
ops.test("NormUnbalance", 1.0e-7, 100, 0)
# ops.test("NormDispIncr", 1.0e-7, 100, 0)

ops.system("UmfPack")


# Creating DOF numberer
ops.numberer("RCM")

# Creating constraints handler
ops.constraints("Plain")

# Creating algorithm
# ops.algorithm("Newton")
# ops.algorithm("Linear")
# ops.algorithm("NewtonLineSearch")
ops.algorithm("NewtonLineSearch", 'type', 'Bisection')


# Create recorder
# ops.recorder('Node', '-file', 'Node12_Z.out', '-closeOnWrite', '-time', '-node', *[12], '-dof', *[3], *['disp'])

# Compute cycle
def computeCycle():
    D0=0
    loadFactor_max=0
    for j in range(1, len(t)):
        delta = uTip[j]-D0
        # print("delta = ", delta)

        # Creating integrator
        ops.integrator("DisplacementControl", 12, 3, delta)  

        # Create analysis type
        ops.analysis("Static")

        if (ops.analyze(1) != 0):
            exit()
        # exit()
        # print("disp12 = ", ops.nodeDisp(12, 3),'\n')
        loadFactor=ops.getLoadFactor(1)
        loadFactor_max=max(loadFactor,loadFactor_max)
        # print("loadFactor = ", loadFactor)
        D0=uTip[j]
    return loadFactor_max

data=np.array([
    [0                  , 106.49681528662421],
    [1186.9436201780336 , 104.39490445859873],
    [2373.887240356067  , 102.05944798301488],
    [3560.830860534101  , 100.19108280254778],
    [9495.548961424298  , 98.32271762208069],
    [16617.210682492558 , 97.38853503184714],
    [27299.703264094976 , 96.4543524416136],
    [46290.80118694363  , 95.75371549893845],
    [68842.72997032641  , 95.28662420382166],
    [89020.7715133531   , 95.05307855626327],
    [108011.86943620178 , 95.05307855626327],
    [127002.96735905044 , 94.5859872611465],
    [148367.95252225522 , 94.5859872611465],
    [172106.824925816   , 94.11889596602973],
    [205341.24629080118 , 93.65180467091295],
    [231454.0059347181  , 93.41825902335458],
    [263501.4836795252  , 93.18471337579618],
    [316913.94658753707 , 92.7176220806794],
    [351335.3115727003  , 92.25053078556265],
    [376261.1275964392  , 91.54989384288749],
    [394065.28189910983 , 91.31634819532908],
    [420178.0415430267  , 90.84925690021232],
    [468842.7299703264  , 89.91507430997876],
    [499703.2640949555  , 89.21443736730362],
    [530563.7982195846  , 88.51380042462847],
    [562611.2759643918  , 87.57961783439491],
    [594658.7537091989  , 86.64543524416136],
    [626706.231454006   , 84.77707006369428],
    [658753.7091988132  , 82.90870488322719],
    [690801.1869436203  , 81.50743099787687],
    [721661.7210682493  , 78.47133757961785],
    [752522.2551928784  , 74.03397027600849],
    [779821.9584569733  , 68.89596602972401],
    [798813.0563798221  , 65.3927813163482]
    ])

cycles_data=data[:,0]
force_data=data[:,1]
steps=np.diff(cycles_data)
midSteps=np.diff(np.linspace(1,1186,50))
steps=np.sort(np.concatenate([[1]*10,[10]*10,[100]*10,[200]*10,[500]*10,[1000]*10,midSteps,steps,[10000]*100]))

# Computing loadFactors
nSteps=len(steps)
loadFactors=np.zeros(nSteps,dtype=np.float64)
cycles=np.zeros(nSteps)
for i in range(nSteps):
    print("Step = ",i, " of ",nSteps)
    loadFactor_max=computeCycle()
    print('loadFactor = ',loadFactor_max)
    loadFactors[i]=loadFactor_max
    nCycles = steps[i]
    ops.setParameter('-val', int(nCycles), "-ele", 1, "advanceDamageState")
    ops.setParameter('-val', 0, "-ele", 1, "resetMaxStress")
    if i==0:
        cycles[i]=nCycles
    else:
        cycles[i]=cycles[i-1]+nCycles

import matplotlib
matplotlib.rc('axes.formatter', useoffset=False)

print(loadFactors)

# plt.rcParams['axes.formatter.useoffset'] = False
plot(cycles,loadFactors,'or')
plot(cycles_data,force_data,'-b')

# ticklabel_format(useOffset=False)
show()
print("Done")



