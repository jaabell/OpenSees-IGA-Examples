#  IGA CANTILEVER PLATE. GEOMETRY OBTAINED FROM NGUYEN'S igaCicleBendingStrip2D


import numpy as np
import opensees as ops
from math import *

# Geomgl utilities for visualization and surface manipulation
from geomdl import NURBS, compatibility, operations, knotvector
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



La = 10;    #
Lb = 1.0;    #
mm = 1.0 / 1000.  # m


ops.wipe()
ops.model('basic', '-ndm', 3, '-ndf', 3)

# These are given in v,u
controlPts = np.array([
[0.00000  , 0.00000 , 0.00000 , 1],
[0.00000  , 1.00000 , 0.00000 , 1],
[0.31250  , 0.00000 , 0.00000 , 1],
[0.31250  , 1.00000 , 0.00000 , 1],
[0.93750  , 0.00000 , 0.00000 , 1],
[0.93750  , 1.00000 , 0.00000 , 1],
[1.87500  , 0.00000 , 0.00000 , 1],
[1.87500  , 1.00000 , 0.00000 , 1],
[3.12500  , 0.00000 , 0.00000 , 1],
[3.12500  , 1.00000 , 0.00000 , 1],
[4.37500  , 0.00000 , 0.00000 , 1],
[4.37500  , 1.00000 , 0.00000 , 1],
[5.62500  , 0.00000 , 0.00000 , 1],
[5.62500  , 1.00000 , 0.00000 , 1],
[6.87500  , 0.00000 , 0.00000 , 1],
[6.87500  , 1.00000 , 0.00000 , 1],
[8.12500  , 0.00000 , 0.00000 , 1],
[8.12500  , 1.00000 , 0.00000 , 1],
[9.06250  , 0.00000 , 0.00000 , 1],
[9.06250  , 1.00000 , 0.00000 , 1],
[9.68750  , 0.00000 , 0.00000 , 1],
[9.68750  , 1.00000 , 0.00000 , 1],
[10.00000 , 0.00000 , 0.0000  , 1],
[10.00000 , 1.00000 , 0.0000  , 1]
])


patchTag = 1
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
surf.knotvector_u = knotvector.generate(surf.degree_u, surf.ctrlpts_size_u)
surf.knotvector_v = knotvector.generate(surf.degree_v, surf.ctrlpts_size_v)


# Visualize surface
surfVisualize(surf, hold=True)


# nDMaterial ElasticIsotropic $nDtag_elastic $elasticidad_probeta
# $poisson_probeta
E1 = 1.2e7  # Young's modulus N/m^2
E2 = E1
nu = 0.3  # Poisson's ratio
rho = 2.0e2  # kg/m^3


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

# matTags = [3, 4, 3, 4, 3]
# thickness = [10. * mm, 10. * mm, 10. * mm, 10. * mm, 10. * mm]
# θ = [0 * deg2rad, 45 * deg2rad, 90 * deg2rad, -45 * deg2rad, 0 * deg2rad]

matTags = [3]
thickness = [10. * mm]
θ = [0 * deg2rad]

I = (Lb * (sum(thickness)**3)) / 12.0

gFact = [0.0, 0.0, 0.0]


Nlayers = len(θ)

controlPts = surf.ctrlpts2d[:]  # Given in v,u
controlPts = np.array(compatibility.flip_ctrlpts2d(controlPts))  # Flipping to u,v


print("controlPts.tolist(): ", controlPts.tolist())

nodeStartTag=1

ops.IGA("Patch", patchTag, nodeStartTag, surf.degree_u, surf.degree_v, surf.ctrlpts_size_u, surf.ctrlpts_size_v,
        "-type", "KLShell",
        # "-nonLinearGeometry", 0,
        "-planeStressMatTags", *matTags,
        "-gFact", *gFact,
        "-theta", *θ,
        "-thickness", *thickness,
        "-uKnot", *surf.knotvector_u, "-vKnot", *surf.knotvector_v, "-controlPts", *controlPts.flatten())

# exit()
nPoints = surf.ctrlpts_size_u * surf.ctrlpts_size_v
print("surf.ctrlpts_size_u: ", surf.ctrlpts_size_u)
print("surf.ctrlpts_size_v: ", surf.ctrlpts_size_v)


fixedNodes = [1, 2, 13, 14]
for n in ops.getNodeTags():
    if n in fixedNodes:
        ops.fix(n, 1, 1, 1)
    # else:
    #     ops.fix(n, 1, 0, 0)


# equalDOFnodes_master = np.arange(2 * surf.ctrlpts_size_u + 1, nPoints, surf.ctrlpts_size_u)
masterNodes = [3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
retainedNodes = [15, 16, 17, 18, 19, 20, 21, 22, 23, 24]

for i in range(len(masterNodes)):
    masterNode = masterNodes[i]
    retainedNode = retainedNodes[i]
    ops.equalDOF(int(masterNode), int(retainedNode), 1)


print("\n\n\nPRINTING DOMAIN-----------------------")
ops.printModel()
print("\n\n\nDONE PRINTING DOMAIN-----------------------")

# exit()


# ------------------------------
# Start of analysis generation
# ------------------------------

# create TimeSeries
# ops.timeSeries("Constant", 1)
ops.timeSeries("Linear", 1)

# create a plain load pattern
ops.pattern("Plain", 1, 1)


print("Loading nodes")

I      = Lb*(sum(thickness)**3)/12.0
G=E1/2/(1+nu)
twistAngle = np.pi*2
It=Lb*(sum(thickness)**3)/3.0

Ip = Lb*sum(thickness)/12*(Lb**2+sum(thickness)**2)

print("Ip: ", Ip)
print("It: ", It)

moment = twistAngle*G*It/La; 
# moment = twistAngle*G*Ip/La; 


nSteps = 250
Pz=moment/Lb

followerLoadsPos = [0.0, 0.0, Pz] 
followerLoadsNeg = [0.0, 0.0, -Pz]

ops.eleLoad("-ele", 1, "-type", "-IGAFollowerLoad", 1.0, 1.0, *followerLoadsPos)
ops.eleLoad("-ele", 1, "-type", "-IGAFollowerLoad", 1.0, 0.0, *followerLoadsNeg)

print("Finished loading nodes")


print("Starting analysis")

# Create test
ops.test("NormDispIncr", 1.0e-3, 300, 1)
# ops.test("NormUnbalance", 1.0e-3, 300, 1)
# ops.test("EnergyIncr", 1.0e-4, 300, 1) # This got to 22.7334 (more than π/2)

# ops.algorithm("Linear")
# ops.algorithm("Newton")
# ops.algorithm("SecantNewton")
# ops.algorithm("NewtonLineSearch",True)
ops.algorithm("NewtonLineSearch",'-type', 'Bisection')
# ops.algorithm("ModifiedNewton")
# ops.algorithm("KrylovNewton")


# create SOE
ops.system("UmfPack")

# create DOF number
# ops.numberer("Plain")
ops.numberer("RCM")

# create constraint handler
ops.constraints("Plain")

# create integrator
ops.integrator("LoadControl", 1.0/nSteps)

# create analysis object
ops.analysis("Static")


# perform the analysis
import matplotlib.pyplot as plt
data = np.zeros((nSteps + 1, 3))

for j in range(nSteps):
  print("=================================")
  print(f"Load step {j}")
  print("=================================")
  result = ops.analyze(1)
  if result != 0:
      break
      exit(-1)
  else:
      # Adding deformation to controlPts
      controlPts = surf.ctrlpts2d[:]
      controlPts = compatibility.flip_ctrlpts2d(controlPts)  # Flipping to u,v

      fDef = 1
      i = 1
      for dim in controlPts:
          for point in dim:
              point[:3] += fDef * np.array(ops.nodeDisp(i))
              i += 1

      # Setting control points for surface
      controlPts = compatibility.flip_ctrlpts2d(controlPts)
      controlPts = (np.array(controlPts).reshape(
          nPoints, 4))
      surf.set_ctrlpts(controlPts.tolist(), surf.ctrlpts_size_u, surf.ctrlpts_size_v)

      # Visualize surface
      # if j==nSteps-1 or j%10==0:
      if j==nSteps-1:
        surfVisualize(surf, hold=True)

      controlPts = surf.ctrlpts2d[:]
      controlPts = compatibility.flip_ctrlpts2d(controlPts)  # Flipping to u,v
      i = 1
      for dim in controlPts:
          for point in dim:
              point[:3] -= fDef * np.array(ops.nodeDisp(i))
              i += 1

      # Setting control points for surface
      controlPts = compatibility.flip_ctrlpts2d(controlPts)
      controlPts = (np.array(controlPts).reshape(
          nPoints, 4))
      surf.set_ctrlpts(controlPts.tolist(), surf.ctrlpts_size_u, surf.ctrlpts_size_v)

      print("ops.nodeDisp(12,3): ", ops.nodeDisp(12,3))
      print("ops.nodeDisp(24,3): ", ops.nodeDisp(24,3))


      # # Obtain A and B from gimmeMCK
      # ops.wipeAnalysis()

      # ops.system('FullGeneral')
      # ops.analysis("Transient")
      # ops.integrator('GimmeMCK', 0.0,0.0,1.0)
      # ops.analyze(1,0.0)
      # N = ops.systemSize() # Has to be done after analyze
  
      # K = ops.printA('-ret') # Or use ops.printA('-file','K.out')
      # K = np.array(K) # Convert the list to an array
      # K.shape = (N,N) # Make the array an NxN matrix

      # print("K: ", K)

      # # Return to normal analysis

      # ops.wipeAnalysis()
      # ops.system("FullGeneral")
      # ops.numberer("Plain")
      # ops.constraints("Plain")
      # ops.integrator("LoadControl", 1.0 / nSteps)
      # ops.algorithm("Newton")
      # ops.test("NormDispIncr", 1.0e-5, 100, 1)
      # ops.analysis("Static")

      print("\nNext load step\n")

surfVisualize(surf, hold=True)
# plt.plot(data[:, 0], data[:, 2], '-or')
# plt.plot(data[:, 1], data[:, 2], '-or')
# plt.xlabel('Horizontal Displacement')
# plt.ylabel('Horizontal Load')
# plt.show()

# print("Done")

# elasticSolution = (Pz * (La**3)) / (3 * E1 * I)

# print("elasticSolution: ", elasticSolution)
# print("data[nSteps,0]: ", data[nSteps, 0])


print("Finished analysis")
