# octave wrapper

from oct2py import octave
from math import *
from matplotlib.pylab import *

from scipy.interpolate import interp1d




octave.addpath(octave.genpath('octaveFiles/'))

a = 54.0;
b = 30.0;
refineLevel = 3
refU=1
refV=1
octave.makeCantileverShell(a,b,refineLevel,refU,refV)
# out = octave.makeCantileverShell()

import scipy.io as sio
shell = sio.loadmat('shell.mat')
import numpy

P = int(shell['p'])
Q = int(shell['q'])
uKnot = shell['uKnot'][0]
vKnot = shell['vKnot'][0]
noPtsX = int(shell['noPtsX'])
noPtsY = int(shell['noPtsY'])
weights = shell['weights']
controlPts = shell['controlPts']




# Geomdl utilities for visualization and surface manipulation
from geomdl import NURBS, compatibility, operations, knotvector, multi
from surfVisualize import *


# Create a BSpline surface instance
surf_L = NURBS.Surface()

# Set surface degrees
surf_L.degree_u = P
surf_L.degree_v = Q


# Setting control points for surface
surf_L.set_ctrlpts(controlPts.tolist(), noPtsY, noPtsX)
controlPts = surf_L.ctrlpts2d[:]
controlPts = compatibility.flip_ctrlpts2d(controlPts)
controlPts=numpy.reshape(controlPts,[noPtsX*noPtsY,4])
surf_L.set_ctrlpts(controlPts.tolist(), noPtsX, noPtsY)




# Set knot vectors
surf_L.knotvector_u = uKnot.tolist()
surf_L.knotvector_v = vKnot.tolist()


noPtsX = surf_L.ctrlpts_size_u
noPtsY = surf_L.ctrlpts_size_v

# Visualize surface
surfVisualize(surf_L, hold=False)




# Creating second surface

a = 50.1;
b = 30.0;
refineLevel = 3
refU=1
refV=1
octave.makeCantileverShell(a,b,refineLevel,refU,refV)
# out = octave.makeCantileverShell()

import scipy.io as sio
shell = sio.loadmat('shell.mat')
import numpy

P = int(shell['p'])
Q = int(shell['q'])
uKnot = shell['uKnot'][0]
vKnot = shell['vKnot'][0]
noPtsX = int(shell['noPtsX'])
noPtsY = int(shell['noPtsY'])
weights = shell['weights']
controlPts = shell['controlPts']




# Create a BSpline surface instance
surf_R = NURBS.Surface()

# Set surface degrees
surf_R.degree_u = P
surf_R.degree_v = Q


# Setting control points for surface
surf_R.set_ctrlpts(controlPts.tolist(), noPtsY, noPtsX)
controlPts = surf_R.ctrlpts2d[:]
controlPts = compatibility.flip_ctrlpts2d(controlPts)
controlPts=numpy.reshape(controlPts,[noPtsX*noPtsY,4])
surf_R.set_ctrlpts(controlPts.tolist(), noPtsX, noPtsY)


# Set knot vectors
surf_R.knotvector_u = uKnot.tolist()
surf_R.knotvector_v = vKnot.tolist()


noPtsX = surf_R.ctrlpts_size_u
noPtsY = surf_R.ctrlpts_size_v


# Translating surface

operations.translate(surf_R,[54,0,0],inplace=True)
# Visualize surface
surfVisualize(surf_R, hold=False)


# Creating bending strip
from edgeHandler import *
interfacePoints, nodesOnLeftPatch = edgeGetter(surf_L, "10", "11")

nodesOnRightPatch = edgeGetter(surf_R, "00", "01")[1]

bendingStrip = makeBendingStrip(
    nodesOnLeftPatch, interfacePoints, nodesOnRightPatch)

# Creating container for multipatches

surfList = [surf_L, surf_R, bendingStrip]

container = multi.SurfaceContainer(surfList)

# Visualize surface
mm=1.0/1000
container.sample_size = 5
for surf in container:
    operations.scale(surf,mm,inplace=True)
    surf.evaluate()
# Visualization configuration
# container.vis = VisVTK.VisSurface(ctrlpts=True, legend=False)
container.vis = VisMPL.VisSurface(ctrlpts=True, legend=False)
# container.vis.ctrlpts_offset=0.1

# Render the surface
evalcolor = ["green", 'blue','red']
container.render(evalcolor=evalcolor)


# # PostProcessing

# import scipy.io as sio
# import h5py

# # file = h5py.File('aspa.mpco', 'r')
# file = h5py.File('/home/rodrigo/FelipeElgueta/OpenSees-IGA-Examples/octaveWrapper/testing/benchmarkDamage.mpco', 'r')

# use_damage = True


# surfList = [surf_L, surf_R]

# times = np.arange(0,2119,40)

# print(len(times))

# stages = [1]
# maxTimes = len(times)
# startTime = 0

# for stage in stages:

#     modelStage = f'MODEL_STAGE[{stage}]'

#     patches = file[modelStage]['MODEL']['ELEMENTS']['1-UnknownMovableObject[0:0]'][:, 0] # These are the subdomains -> Iga Patches
#     elements = file[modelStage]['MODEL']['ELEMENTS']['207-UnknownMovableObject[0:0]'][:, 0] # These are the IGA_KLShell elements

#     nodes = file[modelStage]['MODEL']['ELEMENTS']['207-UnknownMovableObject[0:0]'][:, 1:]

#     print('nodes = \n', nodes)

#     numPatches = 1

#     nodesOnPatches = []
#     elementsOnPatches = []

#     continueSearch = True
#     while continueSearch:

#       for i in range(1,len(nodes)):
#         currentElement = elements[i-1]
#         nextElement = elements[i]

#         print("currentElement = ", currentElement)
#         print("nextElement = ", nextElement)

#         foundPatch = True

#         if abs(nextElement - currentElement) > 1:
#           continueSearch = True
#           foundPatch = True
#           numPatches += 1
#           print('breaking in i = ', i)
#           break_i = i
#           nodesOnPatches.append(nodes[:i,:])
#           elementsOnPatches.append(elements[:i])
#           nodes = nodes[i:,:]
#           elements = elements[i:]
#           print('nodes = ', nodes)
#           print('elements = ', elements)
#           break
#         else:
#           foundPatch = False

#       if not foundPatch:
#         continueSearch = False
#         nodesOnPatches.append(nodes)
#         elementsOnPatches.append(elements)

#       else:
#         pass


#     print('numPatches = ', numPatches)

#     # Getting layup info
#     layupData = file[modelStage]['MODEL']['SECTION_ASSIGNMENTS']
#     layupData_keys = layupData.keys()
#     numLayups = len(layupData)
#     print('numLayups = ', numLayups)
#     print('layupData_keys = ', layupData_keys)

#     zLocs_list = []

#     for key in layupData_keys:
#         print('key = ', key)

#         layupData_here = layupData[key]

#         zLocs_layups = layupData_here['FIBER_DATA'][:,1] # Shape (numMaterial x 3 (algo, zLoc, thickness))

#         print('zLocs_layups = ', zLocs_layups)

#         zLocs_list.append(zLocs_layups)




#     # Getting nodeLimits in order to know which nodes are from which patch
#     nodeLimits = []
#     for nodes in nodesOnPatches:
#       nodeLimits.append([nodes[0,0],nodes[-1,-1]])




#     nTimes = len(file[modelStage]['RESULTS']['ON_NODES']['DISPLACEMENT']['DATA'])

#     nTimes = min(maxTimes,nTimes) # Just so the data is not too much

#     # for time in range(nTimes):
#     for time in range(len(times)):



#       step = f'STEP_{times[time] + startTime}' 

#       print(step)

#       controlPts = file[modelStage]['MODEL']['NODES']['COORDINATES'][:, :] # All the control points / nodes in the model
#       U = file[modelStage]['RESULTS']['ON_NODES']['DISPLACEMENT']['DATA'][step][:,:] # All the displacements in the model

#       nResults = len(file[modelStage]['RESULTS']['ON_ELEMENTS']['section.fiber.stress']) # Doesn't have to be stress, just need the number of files

#       # The files are different because they are not the same size (nLayers)
     
#       dictResults = {}

#       for nFile in range(nResults):
#         # Getting current file tag

#         actualFile = f'207-UnknownMovableObject[0:0:{nFile}]'

#         multiplicity = file[modelStage]['RESULTS']['ON_ELEMENTS']['section.fiber.stress'][actualFile]['META']['MULTIPLICITY'][0,0]


#         dictResults[multiplicity] = actualFile
    

#       controlPts_patches = []
#       sigma_patches = []
#       strain_patches =[]
#       damage_patches = []
#       U_patches=[]
#       Nlayers_patches=[]


#       # Looping through the 'patches'
#       prevNlayers = len(zLocs_list[0])
#       for i in range(len(nodeLimits)):

#         zLocs = zLocs_list[i]

#         Nlayers = len(zLocs)
#         Nlayers_patches.append(Nlayers)

#         # This is [nlayers1, nlayers2, nlayers3,.... etc]
#         sizesFile = dictResults.keys()

#         # for Nlayers_file in sizesFile:
#         #     print('Nlayers_file = ', Nlayers_file)  
#         #     if Nlayers == Nlayers_file

   

#         # Need the specific file on strains and stresses depending on the size
#         actualFile = dictResults[Nlayers]
#         print('actualFile = ', actualFile)


#         sigma_global = file[modelStage]['RESULTS']['ON_ELEMENTS']['section.fiber.stress'][actualFile]['DATA'][step][:,:] # Shape (numElements x Nlayers*noGps*3)
#         strain_global = file[modelStage]['RESULTS']['ON_ELEMENTS']['section.fiber.strain'][actualFile]['DATA'][step][:,:] # Shape (numElements x Nlayers*noGps*3)
        
#         if use_damage:
#             damage_global = file[modelStage]['RESULTS']['ON_ELEMENTS']['section.fiber.damagestate'][actualFile]['DATA'][step][:,:] # Shape (numElements x Nlayers*noGps*9)

#         nodeStart, nodeEnd = nodeLimits[i]
#         rangeNodes = np.arange(nodeStart,nodeEnd+1,1) - 1
#         grid = np.ix_(rangeNodes,range(3))

#         controlPts_here = controlPts[grid]
#         U_here = U[grid]
#         shapeU_here = np.shape(U_here)
#         U_here = np.reshape(U_here,shapeU_here[0]*shapeU_here[1])


#         elementsThisPatch = elementsOnPatches[i]

#         controlPts_patches.append(controlPts_here)
#         U_patches.append(U_here)

#         # print('elementsThisPatch = ', elementsThisPatch)

#         if Nlayers != prevNlayers or i==0:
#             # I'm on a different file
#             newFile = True
#         else:
#             # Same file
#             newFile = False



#         # Have to do a new grid
#         if newFile:
#           # I'm in the first patch
#           numEleThisPatch = len(elementsThisPatch)
#           id_elements = np.arange(0,numEleThisPatch,1)

#         else:
#           elementsPrevPatch = elementsOnPatches[i-1]
#           numEleThisPatch = len(elementsThisPatch)
#           numElePrevPatch = len(elementsPrevPatch)
#           id_elements = np.arange(numElePrevPatch, numElePrevPatch+numEleThisPatch,1)

#         # print('id_elements = ', id_elements)
#         print('len(id_elements) = ', len(id_elements))


#         nData = len(sigma_global[0])
#         grid = np.ix_(id_elements,range(nData))

#         if use_damage:
#             nData_damage = len(damage_global[0])
#             grid_damage = np.ix_(id_elements,range(nData_damage))

#         # continue





#         # Getting stresses and strains from elements in each layer
#         sigma = sigma_global[grid]
#         strain = strain_global[grid]
#         if use_damage:
#             damage = damage_global[grid_damage]

#         # print('sigma = ', sigma)


#         shapeSigma=np.shape(sigma) # numEle x 3*nGauss
#         shapeStrain=np.shape(strain)
#         if use_damage:
#             shapeDamage=np.shape(damage)

#         numEle = shapeSigma[0] 

#         # print('numEle = ', numEle)

#         # Ranges for stress and strain
#         range0 = np.arange(0,shapeSigma[1],3)
#         range1 = range0+1
#         range2 = range0+2

#         # Ranges for damage
#         if use_damage:
#             range0_damage = np.arange(0,shapeDamage[1],9)
#             range1_damage = range0_damage+1
#             range2_damage = range0_damage+2
#             range3_damage = range0_damage+3
#             range4_damage = range0_damage+4
#             range5_damage = range0_damage+5
#             range6_damage = range0_damage+6
#             range7_damage = range0_damage+7
#             range8_damage = range0_damage+8

#         # print('range0 = ', range0)

#         # Getting quantities in each direction for all layers
#         sigmaXX = sigma[:,np.ix_(range0)]
#         sigmaYY = sigma[:,np.ix_(range1)]
#         sigmaXY = sigma[:,np.ix_(range2)]
#         strainXX = strain[:,np.ix_(range0)]
#         strainYY = strain[:,np.ix_(range1)]
#         strainXY = strain[:,np.ix_(range2)]

#         # Getting damage variables (9 in total)
#         if use_damage:
#             damage_0 = damage[:,np.ix_(range0_damage)]
#             damage_1 = damage[:,np.ix_(range1_damage)]
#             damage_2 = damage[:,np.ix_(range2_damage)]
#             damage_3 = damage[:,np.ix_(range3_damage)]
#             damage_4 = damage[:,np.ix_(range4_damage)]
#             damage_5 = damage[:,np.ix_(range5_damage)]
#             damage_6 = damage[:,np.ix_(range6_damage)]
#             damage_7 = damage[:,np.ix_(range7_damage)]
#             damage_8 = damage[:,np.ix_(range8_damage)]


#         # print('sigmaXX = ', sigmaXX)

#         nGauss_3 = int(shapeSigma[1]/Nlayers)
#         if use_damage:
#             nGauss_9 = int(shapeDamage[1]/Nlayers)

#         # print('nGauss_3 = ', nGauss_3)

#         sigma = np.zeros([Nlayers,numEle,nGauss_3])
#         strain = np.zeros([Nlayers,numEle,nGauss_3])
#         if use_damage:
#             damage = np.zeros([Nlayers,numEle,nGauss_9])



#         for iLayer in range(Nlayers):

#           rangeLayer = np.arange(iLayer,int(nGauss_3/3)*Nlayers,Nlayers)
#           if use_damage:
#             rangeLayer_damage = np.arange(iLayer,int(nGauss_9/9)*Nlayers,Nlayers)
#           # if Nlayers == 1:
#           #   rangeLayer = np.arange(iLayer,int(nGauss_3/3),Nlayers)
#           # else:
#           #   rangeLayer = np.arange(iLayer,nGauss_3,Nlayers)

#           # print('rangeLayer = ', rangeLayer)

#           # Getting stresses and strains corresponding to actual layer 
#           sigmaXX_layer = sigmaXX[:,:,np.ix_(rangeLayer)]
#           sigmaYY_layer = sigmaYY[:,:,np.ix_(rangeLayer)]
#           sigmaXY_layer = sigmaXY[:,:,np.ix_(rangeLayer)]


#           # print("sigmaXX_layer = ", sigmaXX_layer)
#           # print("sigmaYY_layer = ", sigmaYY_layer)
#           # print("sigmaXY_layer = ", sigmaXY_layer)
#           # exit()


#           strainXX_layer = strainXX[:,:,np.ix_(rangeLayer)]
#           strainYY_layer = strainYY[:,:,np.ix_(rangeLayer)]
#           strainXY_layer = strainXY[:,:,np.ix_(rangeLayer)]

#           # Getting damage variables corresponding to actual layer
#           if use_damage:
#             damage0_layer = damage_0[:,:,np.ix_(rangeLayer_damage)]
#             damage1_layer = damage_1[:,:,np.ix_(rangeLayer_damage)]
#             damage2_layer = damage_2[:,:,np.ix_(rangeLayer_damage)]
#             damage3_layer = damage_3[:,:,np.ix_(rangeLayer_damage)]
#             damage4_layer = damage_4[:,:,np.ix_(rangeLayer_damage)]
#             damage5_layer = damage_5[:,:,np.ix_(rangeLayer_damage)]
#             damage6_layer = damage_6[:,:,np.ix_(rangeLayer_damage)]
#             damage7_layer = damage_7[:,:,np.ix_(rangeLayer_damage)]
#             damage8_layer = damage_8[:,:,np.ix_(rangeLayer_damage)]




#           shapeSigma = np.shape(sigmaXX_layer)
#           shapeStrain = np.shape(strainXX_layer)
#           if use_damage:
#             shapeDamage = np.shape(damage0_layer)

#           # print('shapeSigma = ', shapeSigma)

#           sigmaXX_layer = np.reshape(sigmaXX_layer,[1,shapeSigma[0],shapeSigma[-1]])
#           sigmaYY_layer = np.reshape(sigmaYY_layer,[1,shapeSigma[0],shapeSigma[-1]])
#           sigmaXY_layer = np.reshape(sigmaXY_layer,[1,shapeSigma[0],shapeSigma[-1]])

#           strainXX_layer = np.reshape(strainXX_layer,[1,shapeStrain[0],shapeStrain[-1]])
#           strainYY_layer = np.reshape(strainYY_layer,[1,shapeStrain[0],shapeStrain[-1]])
#           strainXY_layer = np.reshape(strainXY_layer,[1,shapeStrain[0],shapeStrain[-1]])

#           if use_damage:
#             damage0_layer = np.reshape(damage0_layer,[1,shapeDamage[0],shapeDamage[-1]])
#             damage1_layer = np.reshape(damage1_layer,[1,shapeDamage[0],shapeDamage[-1]])
#             damage2_layer = np.reshape(damage2_layer,[1,shapeDamage[0],shapeDamage[-1]])
#             damage3_layer = np.reshape(damage3_layer,[1,shapeDamage[0],shapeDamage[-1]])
#             damage4_layer = np.reshape(damage4_layer,[1,shapeDamage[0],shapeDamage[-1]])
#             damage5_layer = np.reshape(damage5_layer,[1,shapeDamage[0],shapeDamage[-1]])
#             damage6_layer = np.reshape(damage6_layer,[1,shapeDamage[0],shapeDamage[-1]])
#             damage7_layer = np.reshape(damage7_layer,[1,shapeDamage[0],shapeDamage[-1]])
#             damage8_layer = np.reshape(damage8_layer,[1,shapeDamage[0],shapeDamage[-1]])



#           # range0 = np.arange(0,int(nGauss_3/3),3)
#           # Ranges for stresses and strains
#           range0 = np.arange(0,int(nGauss_3),3)
#           range1 = np.arange(1,int(nGauss_3),3)
#           range2 = np.arange(2,int(nGauss_3),3)

#           # Ranges for damage
#           if use_damage:
#             range0_damage = np.arange(0,int(nGauss_9),9)
#             range1_damage = np.arange(1,int(nGauss_9),9)
#             range2_damage = np.arange(2,int(nGauss_9),9)
#             range3_damage = np.arange(3,int(nGauss_9),9)
#             range4_damage = np.arange(4,int(nGauss_9),9)
#             range5_damage = np.arange(5,int(nGauss_9),9)
#             range6_damage = np.arange(6,int(nGauss_9),9)
#             range7_damage = np.arange(7,int(nGauss_9),9)
#             range8_damage = np.arange(8,int(nGauss_9),9)


#           # Grids for stresses and strains
#           grid0 = np.ix_([iLayer],range(numEle),range0)
#           grid1 = np.ix_([iLayer],range(numEle),range1)
#           grid2 = np.ix_([iLayer],range(numEle),range2)


#           # Grids for damage
#           if use_damage:
#             grid0_damage = np.ix_([iLayer],range(numEle),range0_damage)
#             grid1_damage = np.ix_([iLayer],range(numEle),range1_damage)
#             grid2_damage = np.ix_([iLayer],range(numEle),range2_damage)
#             grid3_damage = np.ix_([iLayer],range(numEle),range3_damage)
#             grid4_damage = np.ix_([iLayer],range(numEle),range4_damage)
#             grid5_damage = np.ix_([iLayer],range(numEle),range5_damage)
#             grid6_damage = np.ix_([iLayer],range(numEle),range6_damage)
#             grid7_damage = np.ix_([iLayer],range(numEle),range7_damage)
#             grid8_damage = np.ix_([iLayer],range(numEle),range8_damage)


#           # Assigning stresses and strains to the arrays
#           sigma[grid0] = sigmaXX_layer
#           sigma[grid1] = sigmaYY_layer
#           sigma[grid2] = sigmaXY_layer
#           strain[grid0] = strainXX_layer
#           strain[grid1] = strainYY_layer
#           strain[grid2] = strainXY_layer

#           # Assigning damage variables to array
#           if use_damage:
#             damage[grid0_damage] = damage0_layer
#             damage[grid1_damage] = damage1_layer
#             damage[grid2_damage] = damage2_layer
#             damage[grid3_damage] = damage3_layer
#             damage[grid4_damage] = damage4_layer
#             damage[grid5_damage] = damage5_layer
#             damage[grid6_damage] = damage6_layer
#             damage[grid7_damage] = damage7_layer
#             damage[grid8_damage] = damage8_layer


#         sigma_patches.append(sigma)
#         strain_patches.append(strain)
#         if use_damage:
#             damage_patches.append(damage)


#         prevNlayers = Nlayers


#       # Until here I have the controlPts for each patch and the U vector for each patch (I hope)

    
#       for (surf,i) in zip(surfList,range(len(surfList))):
#         surf = surfList[i]

#         shell={}

#         Nlayers = Nlayers_patches[i]
#         zLocs = zLocs_list[i]


#         P = surf.degree_u
#         Q = surf.degree_v

#         noPtsX  = surf.ctrlpts_size_u
#         noPtsY  = surf.ctrlpts_size_v


#         uKnot = surf.knotvector_u
#         vKnot = surf.knotvector_v

#         # Recover info for the current patch from the lists
#         controlPts = controlPts_patches[i]
#         weights = [1]*len(controlPts)
#         U = U_patches[i]
#         sigma = sigma_patches[i]
#         strain = strain_patches[i]
#         if use_damage:
#             damage = damage_patches[i]

#         print('np.shape(sigma) = ', np.shape(sigma))


#         shell['p'] = P
#         shell['q'] = Q 
#         shell['uKnot'] = uKnot 
#         shell['vKnot'] = vKnot 
#         shell['noPtsX'] = noPtsX 
#         shell['noPtsY'] = noPtsY 
#         # shell['weights'] = np.transpose(weights)
#         shell['weights'] = weights
#         shell['controlPts'] = controlPts
#         shell['u'] = U 
#         shell['stresses'] = sigma
#         shell['strains'] = strain
#         if use_damage:
#             shell['damage'] = damage
#         shell['nLayers'] = Nlayers
#         shell['zLocs'] = zLocs


#         print('damage = \n', damage,'\n')


#         matName = f'shell_{i}_{time}.mat'
#         sio.savemat(matName,shell)

#         print("Done saving .mat = ", matName)










import opensees as ops

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

vonPaepeParams = [1e2*E1, 1e2*E2, nu12, nu21, 1e2*G12, rho, Xt, Xc,
                  Yt, Yc, S, c1, c2, c3, c4, c5, c6, c7, c8, c9, b]
tagPlaneStress2 = 2
ops.nDMaterial("VonPapaDamage", tagPlaneStress2, *vonPaepeParams)

tagNDmat2 = 3
ops.nDMaterial("ElasticIsotropic", tagNDmat2, 1e3*E1, 0, rho)
tagPlaneStress3 = 4
ops.nDMaterial("PlaneStress", tagPlaneStress3, tagNDmat2)

materialTags = [tagPlaneStress1, tagPlaneStress2, tagPlaneStress3]


# Laminated parameters

deg2rad = pi / 180

nLayers = 8
totalThickness = 2.72*mm
thick_i = totalThickness/nLayers
thickness = [thick_i]*nLayers
θ = deg2rad*np.array([0]*nLayers)
matTags = [materialTags[0]]*nLayers


names = ["surf_L", "surf_La", 'bendingStrip']
shellType = 'KLShell'
patchTag = 1
nodeStartTag = 1

nodesMap = []

for i in range(len(container)):
    surf = container[i]
    name = names[i]

    print("Creating patch ", name)

    if name == 'surf_La':
        # matTags = [materialTags[1]]
        # θ = [0]
        # thickness = [totalThickness]
        pass
    elif name == 'bendingStrip':
        shellType = "KLShell_BendingStrip"
        matTags = [materialTags[2]]
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
    ops.IGA("SurfacePatch", patchTag, nodeStartTag, surf.degree_u, surf.degree_v, surf.ctrlpts_size_u, surf.ctrlpts_size_v,
            "-type", shellType,
            # "-nonLinearGeometry", 0,
            "-planeStressMatTags", *matTags,
            "-theta", *θ,
            "-thickness", *thickness,
            "-uKnot", *surf.knotvector_u, "-vKnot", *surf.knotvector_v, "-controlPts", *controlPts.flatten())

    # Get the nodes on current patch
    nodesMap.append(
        np.arange(nodeStartTag, ops.getNodeTags()[-1] + 1).tolist())

    # Update patchTag, nodeStartTag
    lastElTag = ops.getEleTags()[-1]
    lastNodeTag = ops.getNodeTags()[-1]
    patchTag = lastElTag + 1
    nodeStartTag = lastNodeTag + 1



print("\n\n\nPRINTING DOMAIN-----------------------")
ops.printModel()
print("\n\n\nDONE PRINTING DOMAIN-----------------------")


# Getting surfR info
print("\nGetting surf_R info\n")
nPoints_R = surf_R.ctrlpts_size_u * surf_R.ctrlpts_size_v
middleNodes_R = np.arange(1, nPoints_R + 1, surf_R.ctrlpts_size_u) + nPoints_R 
# middleNodes_R = np.arange(1, nPoints_R + 1, surf_R.ctrlpts_size_u)
nextToMiddleNodes_R = middleNodes_R + 1

# Getting surfL info
print("\nGetting surf_L info\n")
nPoints_L = surf_L.ctrlpts_size_u * surf_L.ctrlpts_size_v
middleNodes_L = np.arange(surf_L.ctrlpts_size_u, nPoints_L+1, surf_L.ctrlpts_size_u)
# middleNodes_L = np.arange(surf_L.ctrlpts_size_u, nPoints_L+1, surf_L.ctrlpts_size_u) 
nextToMiddleNodes_L = middleNodes_L - 1

# EqualDofing overlapping points (constrainer Right, retained Left)

for i in range(len(middleNodes_R)):
    retainedNode = int(middleNodes_R[i])
    constrainedNode = int(middleNodes_L[i])
    ops.equalDOF(retainedNode, constrainedNode, 1, 2, 3)




# EqualDofing bending strip with constrainerNodes

# Retaining Nodes
firstRow_r = nextToMiddleNodes_R
secondRow_r = middleNodes_R
thirdRow_r = nextToMiddleNodes_L



# Constrained nodes
nodesOnBendingStrip = nodesMap[2]


firstRow_c = []
secondRow_c = []
thirdRow_c = []

# firstRow_c=(nodesOnBendingStrip[:len(firstRow_r)+1])
# secondRow_c=(nodesOnBendingStrip[len(firstRow_r):2*len(firstRow_r)])
# thirdRow_c=(nodesOnBendingStrip[2*len(firstRow_r):])

for i in range(int(len(nodesOnBendingStrip)/3)):
    thirdRow_c.append(nodesOnBendingStrip[3 * i ])
    secondRow_c.append(nodesOnBendingStrip[3 * i + 1])
    firstRow_c.append(nodesOnBendingStrip[3 * i + 2])


for i in range(len(firstRow_c)):
    retainedNode = int(firstRow_r[i])
    constrainedNode = int(firstRow_c[i])
    ops.equalDOF(retainedNode, constrainedNode, 1, 2, 3)


for i in range(len(secondRow_c)):
    retainedNode = int(secondRow_r[i])
    constrainedNode = int(secondRow_c[i])
    ops.equalDOF(retainedNode, constrainedNode, 1, 2, 3)



for i in range(len(thirdRow_c)):
    retainedNode = int(thirdRow_r[i])
    constrainedNode = int(thirdRow_c[i])
    ops.equalDOF(retainedNode, constrainedNode, 1, 2, 3)


fixedNodes_L = np.arange(1, nPoints_L + 1, surf_L.ctrlpts_size_u) 
nextToFixedNodes_L = fixedNodes_L + 1

# Creating constraints
for n in ops.getNodeTags():
    n = int(n)
    if n in fixedNodes_L or n in nextToFixedNodes_L:
        ops.fix(n, 1, 1, 1)



ω = 2*pi  # rad/s
tMax = 1  # 5 seconds
deltaT = 2.5e-2
t = np.arange(0, tMax+deltaT, deltaT)

uMax = 30.4*mm/2
# uMax = 1.2*mm/2
uTip = np.sin(ω*t-pi/2)*uMax+uMax

plot(t, 1000*uTip, '-o')
show()

# ------------------------------
# Start of analysis generation
# ------------------------------

# Create MPCO recorder
ops.recorder("mpco","benchmarkDamage",
    "-N","displacement",
    "-E","section.fiber.stress",
    "-E","section.fiber.strain",
    "-E","section.fiber.damagestate",
    )


# create TimeSeries
ops.timeSeries("Linear", 1)

# create a plain load pattern
ops.pattern("Plain", 1, 1)


print("Loading nodes")

nodesToLoad = middleNodes_R + surf_R.ctrlpts_size_u - 1
Pz = 1.0/(2*(len(nodesToLoad)-1))

for n in nodesToLoad:
  n=int(n)
  if n == nodesToLoad[0] or n == nodesToLoad[-1]:
    ops.load(n, 0, 0, Pz)
  else:
    ops.load(n, 0, 0, 2*Pz)



print("Finished loading nodes")


print("Starting analysis")

# Create test
# ops.test("EnergyIncr", 1.0e-5, 100, 0)
ops.test("NormUnbalance", 1.0e-5, 100, 0)
# ops.test("NormDispIncr", 1.0e-5, 100, 0)

# create SOE
ops.system("UmfPack")
# ops.system("FullGeneral")

# create DOF number
ops.numberer("RCM")

# create constraint handler
ops.constraints("Plain")
# ops.constraints("Penalty",1,1)


# ops.algorithm("Linear")
# ops.algorithm("Newton")
# ops.algorithm("SecantNewton")
# ops.algorithm("NewtonLineSearch")
# ops.algorithm("ModifiedNewton")
# ops.algorithm("KrylovNewton")
# ops.algorithm("BFGS")
# ops.algorithm("Broyden")

# # Create analysis type
# ops.analysis("Static")





# Compute cycle

# Adding deformation to control points
fdef = 1e2
for i in range(len(container)):
    surf=container[i]
    nodes=nodesMap[i]
    controlPts = np.array(compatibility.flip_ctrlpts2d(surf.ctrlpts2d[:])).reshape([surf.ctrlpts_size_v * surf.ctrlpts_size_u, 4]).tolist()
    for n in nodes:
        # Get node position in list
        indexN = nodes.index(n)
        point = controlPts[indexN]

        # Add deformation scaled by fdef
        weight=point[3]
        point[0] = (ops.nodeCoord(n)[0] + fdef * ops.nodeDisp(n)[0])*weight/mm
        point[1] = (ops.nodeCoord(n)[1] + fdef * ops.nodeDisp(n)[1])*weight/mm
        point[2] = (ops.nodeCoord(n)[2] + fdef * ops.nodeDisp(n)[2])*weight/mm


    nPoints=surf.ctrlpts_size_u*surf.ctrlpts_size_v
    shape=np.array(compatibility.flip_ctrlpts2d(surf.ctrlpts2d[:])).shape
    controlPts=np.array(controlPts).reshape(shape)
    controlPts=np.array(compatibility.flip_ctrlpts2d(controlPts))

    surf.set_ctrlpts(controlPts.reshape(nPoints,4).tolist(), surf.ctrlpts_size_u, surf.ctrlpts_size_v)



container.sample_size = 10
for surf in container:
    surf.evaluate()
    
# Visualization configuration
# container.vis = VisMPL.VisSurface(ctrlpts=True, legend=False, animate=True)
container.vis = VisVTK.VisSurface(ctrlpts=True, legend=False, animate=True)
# container.vis = VisMPL.VisSurfWireframe(ctrlpts=True, legend=False)
# container.vis.ctrlpts_offset=0.1

# Render the surface
evalcolor = ["green", "green", 'black']
cpcolor=["green","green","black"]
container.render(evalcolor=evalcolor, cpcolor=cpcolor)

# exit()


def computeCycle():
    D0 = 0
    loadFactor_max = 0
    for j in range(1, len(t)):
        delta = uTip[j]-D0
        # print("delta = ", delta)

        # Creating integrator
        ops.integrator("DisplacementControl", int(nodesToLoad[-1]), 3, delta)

        # Creating algorithm
        ops.algorithm("Newton")
        # ops.algorithm("NewtonLineSearch")

        # Create analysis type
        ops.analysis("Static")


        if (ops.analyze(1) != 0):
            print("Analysis failed")
            return 0
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
# steps = np.sort(np.concatenate(
    # [[1]*10, [10]*10, [100]*10, [200]*10, [500]*10, [1000]*10, midSteps, steps, [1000]*10]))
# steps = np.concatenate([[1]*10, [10]*10, [100]*10, [200]*30,[500]*30,[1000]*780])#, [500]*10, [1000]*10, [2000]*10, [5000]*20])
steps = np.concatenate([[1]*10, [10]*10, [100]*10, [200]*30,[500]*30,[1000]*30,[5000]*30,[10000]*40,[5000]*30,[1000]*50])# This works well
steps = np.concatenate([[1]*1, [10]*5, [100]*5, [1000]*10,[20000]*10,[60000]*10,[10000]*30])#, [500]*10, [1000]*10, [2000]*10, [5000]*20])
steps = np.concatenate([[1]*1, [10]*5, [100]*5, [1000]*10,[20000]*10,[30000]*14,[15000]*20])# Funciona bien



print(sum(steps)/798813)
print(sum(steps)/1e5)
print(len(steps))
print("lastElTag = ",lastElTag)


surf= container[0]
p=surf.degree_u
q=surf.degree_v
nGauss=(p+1)*(q+1)
elementsWithDamage=np.arange(2,65+1,dtype=int)
nElems=len(elementsWithDamage)
nCapas=nLayers

NJUMPVECTOR=np.zeros([nElems,nGauss,nCapas,3],dtype=int)

# exit()
# Computing loadFactors
nSteps = len(steps)
# loadFactors = np.zeros(nSteps, dtype=np.float64)
loadFactors = []
# cycles = np.zeros(nSteps)
cycles = []

cycle=1
i=0
while (cycle < 8e5):
# for i in range(nSteps):
    print("\n\nStep = ", i)
    loadFactor_max = computeCycle()
    if loadFactor_max == 0:
        print("cycles = ", cycles)
        print("loadFactors = ", loadFactors)
        break
    print('loadFactor = ', loadFactor_max, '\n')
    # loadFactors[i] = loadFactor_max
    loadFactors.append(loadFactor_max)
    # nCycles = steps[i]
    # nCycles.append(cycle)

    for eleNumber in elementsWithDamage:
        for gp in range(nGauss):
            for capa in range(nCapas):
                eleNumber = int(eleNumber)
                gp=int(gp)
                capa=int(capa)
                NJUMPVECTOR[eleNumber-elementsWithDamage[0], gp, capa,:] = ops.eleResponse(eleNumber, "material_and_layer_number", str(gp), str(capa), 'NJUMP', '0.01')
                # print('NJUMPVECTOR[eleNumber-elementsWithDamage[0], gp, capa, :] = ', NJUMPVECTOR[eleNumber-elementsWithDamage[0], gp, capa, :])
    
    # kwargs = dict(alpha=0.5, bins=100, density=True, stacked=True)
    # hist(NJUMPVECTOR.flatten(),**kwargs)

    hist, bins = np.histogram(NJUMPVECTOR.flatten(), bins=int(1e6), density=True)
    # hist = hist/len(NJUMPVECTOR.flatten())
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2


    cumulativeFreq = np.cumsum(np.diff(bins)*hist) # This is the integral of the distribution
    cumulativeFreq = np.insert(cumulativeFreq,0,[0.0])

    # np.insert(cumulativeFreq,-1,[1.0])
    # np.insert(bins,0,[0.0])
    # cumulativeFreq_max = cumulativeFreq[-1] # This is the maximum of the previous integral
    # cumulativeFreq /= cumulativeFreq_max # This is the normalized integral
    
    # cumulativeFreq /= len(NJUMPVECTOR.flatten()) # This is the normalized integral


    NJUMP_perc = 1e-2
    index_njump_perc = (np.abs(cumulativeFreq-NJUMP_perc)).argmin()
    print('index_njump_perc = ', index_njump_perc)
    print("int(bins[index_njump_perc]) = ", int(bins[index_njump_perc]))

    for j in range(len(cumulativeFreq)):
        print("cumulativeFreq[j] = ", cumulativeFreq[j])
        if cumulativeFreq[j] >= NJUMP_perc:
            index_njump_perc = j-1
            print("int(bins[index_njump_perc]) = ", int(bins[index_njump_perc]))
            break

    print('bins = ',bins)

    nCycles = int(bins[index_njump_perc])

    frequency_inter=interp1d(cumulativeFreq,bins)
    nCycles = int(frequency_inter(NJUMP_perc))
    print('frequency_inter(NJUMP_perc) = ',frequency_inter(NJUMP_perc))

    if nCycles == 0 :
        print("nCycles = 0 ")
        nCycles = 5
        # break


    # print("cumulativeFreq[-1] = ", cumulativeFreq[-1])
    print('nJump = ', nCycles)

    # bar(center, hist, align='center', width=width)

    plot(bins,cumulativeFreq)
    # plot(frequency_inter(cumulativeFreq),cumulativeFreq, '--')
    # plot(nCycles,NJUMP_perc,'or')
    plot(nCycles,NJUMP_perc,'ob')

    # show()
    savefig("{}.png".format(i))
    clf()

   

    # exit()




    # nCycles = center[1]

    # PARAPIPE

    # exit()

    print("Advancing damage state in ", nCycles, "cycles")
    ops.setParameter('-val', int(nCycles), "-ele", 1, "advanceDamageState") # Only calling to patch 1

    # Step integration

    # if nCycles>=100:
    #     for j in range(100):
    #         ops.setParameter('-val', int(nCycles/100), "-ele", 1, "advanceDamageState")
    # else:
    #     ops.setParameter('-val', int(nCycles), "-ele", 1, "advanceDamageState")

    ops.setParameter('-val', 0, "-ele", 1, "resetMaxStress")

    if i == 0:
        # cycles[i] = nCycles
        cycles.append(nCycles)
    else:
        cycles.append(cycles[i-1]+nCycles)

    print("Cycles = ", cycles[i])
    i+=1
    cycle+=nCycles

  

    matplotlib.rc('axes.formatter', useoffset=False)

    print("loadFactors = \n", loadFactors)
    print("Cycles = \n",cycles)

    # plt.rcParams['axes.formatter.useoffset'] = False
    plot(cycles, loadFactors, 'or')
    plot(cycles_data, force_data, '-b')
    savefig("p{}.png".format(i))
    clf()

matplotlib.rc('axes.formatter', useoffset=False)

print("loadFactors = \n", loadFactors)
print("Cycles = \n",cycles)

# plt.rcParams['axes.formatter.useoffset'] = False
plot(cycles, loadFactors, 'or')
plot(cycles_data, force_data, '-b')

print(loadFactors)



# ticklabel_format(useOffset=False)
show()
print("Done")


# PostProcessing

import scipy.io as sio
import h5py

# file = h5py.File('aspa.mpco', 'r')
file = h5py.File('/home/rodrigo/FelipeElgueta/OpenSees-IGA-Examples/octaveWrapper/testing/benchmarkDamage.mpco', 'r')

use_damage = True


surfList = [surf_L, surf_R]

times = np.arange(0,2119,40)

print(len(times))

stages = [1]
maxTimes = len(times)
startTime = 0

for stage in stages:

    modelStage = f'MODEL_STAGE[{stage}]'

    patches = file[modelStage]['MODEL']['ELEMENTS']['1-UnknownMovableObject[0:0]'][:, 0] # These are the subdomains -> Iga Patches
    elements = file[modelStage]['MODEL']['ELEMENTS']['207-UnknownMovableObject[0:0]'][:, 0] # These are the IGA_KLShell elements

    nodes = file[modelStage]['MODEL']['ELEMENTS']['207-UnknownMovableObject[0:0]'][:, 1:]

    print('nodes = \n', nodes)

    numPatches = 1

    nodesOnPatches = []
    elementsOnPatches = []

    continueSearch = True
    while continueSearch:

      for i in range(1,len(nodes)):
        currentElement = elements[i-1]
        nextElement = elements[i]

        print("currentElement = ", currentElement)
        print("nextElement = ", nextElement)

        foundPatch = True

        if abs(nextElement - currentElement) > 1:
          continueSearch = True
          foundPatch = True
          numPatches += 1
          print('breaking in i = ', i)
          break_i = i
          nodesOnPatches.append(nodes[:i,:])
          elementsOnPatches.append(elements[:i])
          nodes = nodes[i:,:]
          elements = elements[i:]
          print('nodes = ', nodes)
          print('elements = ', elements)
          break
        else:
          foundPatch = False

      if not foundPatch:
        continueSearch = False
        nodesOnPatches.append(nodes)
        elementsOnPatches.append(elements)

      else:
        pass


    print('numPatches = ', numPatches)

    # Getting layup info
    layupData = file[modelStage]['MODEL']['SECTION_ASSIGNMENTS']
    layupData_keys = layupData.keys()
    numLayups = len(layupData)
    print('numLayups = ', numLayups)
    print('layupData_keys = ', layupData_keys)

    zLocs_list = []

    for key in layupData_keys:
        print('key = ', key)

        layupData_here = layupData[key]

        zLocs_layups = layupData_here['FIBER_DATA'][:,1] # Shape (numMaterial x 3 (algo, zLoc, thickness))

        print('zLocs_layups = ', zLocs_layups)

        zLocs_list.append(zLocs_layups)




    # Getting nodeLimits in order to know which nodes are from which patch
    nodeLimits = []
    for nodes in nodesOnPatches:
      nodeLimits.append([nodes[0,0],nodes[-1,-1]])




    nTimes = len(file[modelStage]['RESULTS']['ON_NODES']['DISPLACEMENT']['DATA'])

    nTimes = min(maxTimes,nTimes) # Just so the data is not too much

    # for time in range(nTimes):
    for time in times:

      step = f'STEP_{time + startTime}' 

      controlPts = file[modelStage]['MODEL']['NODES']['COORDINATES'][:, :] # All the control points / nodes in the model
      U = file[modelStage]['RESULTS']['ON_NODES']['DISPLACEMENT']['DATA'][step][:,:] # All the displacements in the model

      nResults = len(file[modelStage]['RESULTS']['ON_ELEMENTS']['section.fiber.stress']) # Doesn't have to be stress, just need the number of files

      # The files are different because they are not the same size (nLayers)
     
      dictResults = {}

      for nFile in range(nResults):
        # Getting current file tag

        actualFile = f'207-UnknownMovableObject[0:0:{nFile}]'

        multiplicity = file[modelStage]['RESULTS']['ON_ELEMENTS']['section.fiber.stress'][actualFile]['META']['MULTIPLICITY'][0,0]


        dictResults[multiplicity] = actualFile
    

      controlPts_patches = []
      sigma_patches = []
      strain_patches =[]
      damage_patches = []
      U_patches=[]
      Nlayers_patches=[]


      # Looping through the 'patches'
      prevNlayers = len(zLocs_list[0])
      for i in range(len(nodeLimits)):

        zLocs = zLocs_list[i]

        Nlayers = len(zLocs)
        Nlayers_patches.append(Nlayers)

        # This is [nlayers1, nlayers2, nlayers3,.... etc]
        sizesFile = dictResults.keys()

        # for Nlayers_file in sizesFile:
        #     print('Nlayers_file = ', Nlayers_file)  
        #     if Nlayers == Nlayers_file

   

        # Need the specific file on strains and stresses depending on the size
        actualFile = dictResults[Nlayers]
        print('actualFile = ', actualFile)


        sigma_global = file[modelStage]['RESULTS']['ON_ELEMENTS']['section.fiber.stress'][actualFile]['DATA'][step][:,:] # Shape (numElements x Nlayers*noGps*3)
        strain_global = file[modelStage]['RESULTS']['ON_ELEMENTS']['section.fiber.strain'][actualFile]['DATA'][step][:,:] # Shape (numElements x Nlayers*noGps*3)
        
        if use_damage:
            damage_global = file[modelStage]['RESULTS']['ON_ELEMENTS']['section.fiber.damagestate'][actualFile]['DATA'][step][:,:] # Shape (numElements x Nlayers*noGps*9)

        nodeStart, nodeEnd = nodeLimits[i]
        rangeNodes = np.arange(nodeStart,nodeEnd+1,1) - 1
        grid = np.ix_(rangeNodes,range(3))

        controlPts_here = controlPts[grid]
        U_here = U[grid]
        shapeU_here = np.shape(U_here)
        U_here = np.reshape(U_here,shapeU_here[0]*shapeU_here[1])


        elementsThisPatch = elementsOnPatches[i]

        controlPts_patches.append(controlPts_here)
        U_patches.append(U_here)

        # print('elementsThisPatch = ', elementsThisPatch)

        if Nlayers != prevNlayers or i==0:
            # I'm on a different file
            newFile = True
        else:
            # Same file
            newFile = False



        # Have to do a new grid
        if newFile:
          # I'm in the first patch
          numEleThisPatch = len(elementsThisPatch)
          id_elements = np.arange(0,numEleThisPatch,1)

        else:
          elementsPrevPatch = elementsOnPatches[i-1]
          numEleThisPatch = len(elementsThisPatch)
          numElePrevPatch = len(elementsPrevPatch)
          id_elements = np.arange(numElePrevPatch, numElePrevPatch+numEleThisPatch,1)

        # print('id_elements = ', id_elements)
        print('len(id_elements) = ', len(id_elements))


        nData = len(sigma_global[0])
        grid = np.ix_(id_elements,range(nData))

        if use_damage:
            nData_damage = len(damage_global[0])
            grid_damage = np.ix_(id_elements,range(nData_damage))

        # continue





        # Getting stresses and strains from elements in each layer
        sigma = sigma_global[grid]
        strain = strain_global[grid]
        if use_damage:
            damage = damage_global[grid_damage]

        # print('sigma = ', sigma)


        shapeSigma=np.shape(sigma) # numEle x 3*nGauss
        shapeStrain=np.shape(strain)
        if use_damage:
            shapeDamage=np.shape(damage)

        numEle = shapeSigma[0] 

        # print('numEle = ', numEle)

        # Ranges for stress and strain
        range0 = np.arange(0,shapeSigma[1],3)
        range1 = range0+1
        range2 = range0+2

        # Ranges for damage
        if use_damage:
            range0_damage = np.arange(0,shapeDamage[1],9)
            range1_damage = range0_damage+1
            range2_damage = range0_damage+2
            range3_damage = range0_damage+3
            range4_damage = range0_damage+4
            range5_damage = range0_damage+5
            range6_damage = range0_damage+6
            range7_damage = range0_damage+7
            range8_damage = range0_damage+8

        # print('range0 = ', range0)

        # Getting quantities in each direction for all layers
        sigmaXX = sigma[:,np.ix_(range0)]
        sigmaYY = sigma[:,np.ix_(range1)]
        sigmaXY = sigma[:,np.ix_(range2)]
        strainXX = strain[:,np.ix_(range0)]
        strainYY = strain[:,np.ix_(range1)]
        strainXY = strain[:,np.ix_(range2)]

        # Getting damage variables (9 in total)
        if use_damage:
            damage_0 = damage[:,np.ix_(range0_damage)]
            damage_1 = damage[:,np.ix_(range1_damage)]
            damage_2 = damage[:,np.ix_(range2_damage)]
            damage_3 = damage[:,np.ix_(range3_damage)]
            damage_4 = damage[:,np.ix_(range4_damage)]
            damage_5 = damage[:,np.ix_(range5_damage)]
            damage_6 = damage[:,np.ix_(range6_damage)]
            damage_7 = damage[:,np.ix_(range7_damage)]
            damage_8 = damage[:,np.ix_(range8_damage)]


        # print('sigmaXX = ', sigmaXX)

        nGauss_3 = int(shapeSigma[1]/Nlayers)
        if use_damage:
            nGauss_9 = int(shapeDamage[1]/Nlayers)

        # print('nGauss_3 = ', nGauss_3)

        sigma = np.zeros([Nlayers,numEle,nGauss_3])
        strain = np.zeros([Nlayers,numEle,nGauss_3])
        if use_damage:
            damage = np.zeros([Nlayers,numEle,nGauss_9])



        for iLayer in range(Nlayers):

          rangeLayer = np.arange(iLayer,int(nGauss_3/3)*Nlayers,Nlayers)
          if use_damage:
            rangeLayer_damage = np.arange(iLayer,int(nGauss_9/9)*Nlayers,Nlayers)
          # if Nlayers == 1:
          #   rangeLayer = np.arange(iLayer,int(nGauss_3/3),Nlayers)
          # else:
          #   rangeLayer = np.arange(iLayer,nGauss_3,Nlayers)

          # print('rangeLayer = ', rangeLayer)

          # Getting stresses and strains corresponding to actual layer 
          sigmaXX_layer = sigmaXX[:,:,np.ix_(rangeLayer)]
          sigmaYY_layer = sigmaYY[:,:,np.ix_(rangeLayer)]
          sigmaXY_layer = sigmaXY[:,:,np.ix_(rangeLayer)]


          # print("sigmaXX_layer = ", sigmaXX_layer)
          # print("sigmaYY_layer = ", sigmaYY_layer)
          # print("sigmaXY_layer = ", sigmaXY_layer)
          # exit()


          strainXX_layer = strainXX[:,:,np.ix_(rangeLayer)]
          strainYY_layer = strainYY[:,:,np.ix_(rangeLayer)]
          strainXY_layer = strainXY[:,:,np.ix_(rangeLayer)]

          # Getting damage variables corresponding to actual layer
          if use_damage:
            damage0_layer = damage_0[:,:,np.ix_(rangeLayer_damage)]
            damage1_layer = damage_1[:,:,np.ix_(rangeLayer_damage)]
            damage2_layer = damage_2[:,:,np.ix_(rangeLayer_damage)]
            damage3_layer = damage_3[:,:,np.ix_(rangeLayer_damage)]
            damage4_layer = damage_4[:,:,np.ix_(rangeLayer_damage)]
            damage5_layer = damage_5[:,:,np.ix_(rangeLayer_damage)]
            damage6_layer = damage_6[:,:,np.ix_(rangeLayer_damage)]
            damage7_layer = damage_7[:,:,np.ix_(rangeLayer_damage)]
            damage8_layer = damage_8[:,:,np.ix_(rangeLayer_damage)]




          shapeSigma = np.shape(sigmaXX_layer)
          shapeStrain = np.shape(strainXX_layer)
          if use_damage:
            shapeDamage = np.shape(damage0_layer)

          # print('shapeSigma = ', shapeSigma)

          sigmaXX_layer = np.reshape(sigmaXX_layer,[1,shapeSigma[0],shapeSigma[-1]])
          sigmaYY_layer = np.reshape(sigmaYY_layer,[1,shapeSigma[0],shapeSigma[-1]])
          sigmaXY_layer = np.reshape(sigmaXY_layer,[1,shapeSigma[0],shapeSigma[-1]])

          strainXX_layer = np.reshape(strainXX_layer,[1,shapeStrain[0],shapeStrain[-1]])
          strainYY_layer = np.reshape(strainYY_layer,[1,shapeStrain[0],shapeStrain[-1]])
          strainXY_layer = np.reshape(strainXY_layer,[1,shapeStrain[0],shapeStrain[-1]])

          if use_damage:
            damage0_layer = np.reshape(damage0_layer,[1,shapeDamage[0],shapeDamage[-1]])
            damage1_layer = np.reshape(damage1_layer,[1,shapeDamage[0],shapeDamage[-1]])
            damage2_layer = np.reshape(damage2_layer,[1,shapeDamage[0],shapeDamage[-1]])
            damage3_layer = np.reshape(damage3_layer,[1,shapeDamage[0],shapeDamage[-1]])
            damage4_layer = np.reshape(damage4_layer,[1,shapeDamage[0],shapeDamage[-1]])
            damage5_layer = np.reshape(damage5_layer,[1,shapeDamage[0],shapeDamage[-1]])
            damage6_layer = np.reshape(damage6_layer,[1,shapeDamage[0],shapeDamage[-1]])
            damage7_layer = np.reshape(damage7_layer,[1,shapeDamage[0],shapeDamage[-1]])
            damage8_layer = np.reshape(damage8_layer,[1,shapeDamage[0],shapeDamage[-1]])



          # range0 = np.arange(0,int(nGauss_3/3),3)
          # Ranges for stresses and strains
          range0 = np.arange(0,int(nGauss_3),3)
          range1 = np.arange(1,int(nGauss_3),3)
          range2 = np.arange(2,int(nGauss_3),3)

          # Ranges for damage
          if use_damage:
            range0_damage = np.arange(0,int(nGauss_9),9)
            range1_damage = np.arange(1,int(nGauss_9),9)
            range2_damage = np.arange(2,int(nGauss_9),9)
            range3_damage = np.arange(3,int(nGauss_9),9)
            range4_damage = np.arange(4,int(nGauss_9),9)
            range5_damage = np.arange(5,int(nGauss_9),9)
            range6_damage = np.arange(6,int(nGauss_9),9)
            range7_damage = np.arange(7,int(nGauss_9),9)
            range8_damage = np.arange(8,int(nGauss_9),9)


          # Grids for stresses and strains
          grid0 = np.ix_([iLayer],range(numEle),range0)
          grid1 = np.ix_([iLayer],range(numEle),range1)
          grid2 = np.ix_([iLayer],range(numEle),range2)


          # Grids for damage
          if use_damage:
            grid0_damage = np.ix_([iLayer],range(numEle),range0_damage)
            grid1_damage = np.ix_([iLayer],range(numEle),range1_damage)
            grid2_damage = np.ix_([iLayer],range(numEle),range2_damage)
            grid3_damage = np.ix_([iLayer],range(numEle),range3_damage)
            grid4_damage = np.ix_([iLayer],range(numEle),range4_damage)
            grid5_damage = np.ix_([iLayer],range(numEle),range5_damage)
            grid6_damage = np.ix_([iLayer],range(numEle),range6_damage)
            grid7_damage = np.ix_([iLayer],range(numEle),range7_damage)
            grid8_damage = np.ix_([iLayer],range(numEle),range8_damage)


          # Assigning stresses and strains to the arrays
          sigma[grid0] = sigmaXX_layer
          sigma[grid1] = sigmaYY_layer
          sigma[grid2] = sigmaXY_layer
          strain[grid0] = strainXX_layer
          strain[grid1] = strainYY_layer
          strain[grid2] = strainXY_layer

          # Assigning damage variables to array
          if use_damage:
            damage[grid0_damage] = damage0_layer
            damage[grid1_damage] = damage1_layer
            damage[grid2_damage] = damage2_layer
            damage[grid3_damage] = damage3_layer
            damage[grid4_damage] = damage4_layer
            damage[grid5_damage] = damage5_layer
            damage[grid6_damage] = damage6_layer
            damage[grid7_damage] = damage7_layer
            damage[grid8_damage] = damage8_layer


        sigma_patches.append(sigma)
        strain_patches.append(strain)
        if use_damage:
            damage_patches.append(damage)


        prevNlayers = Nlayers


      # Until here I have the controlPts for each patch and the U vector for each patch (I hope)

    
      for (surf,i) in zip(surfList,range(len(surfList))):
        surf = surfList[i]

        shell={}

        Nlayers = Nlayers_patches[i]
        zLocs = zLocs_list[i]


        P = surf.degree_u
        Q = surf.degree_v

        noPtsX  = surf.ctrlpts_size_u
        noPtsY  = surf.ctrlpts_size_v


        uKnot = surf.knotvector_u
        vKnot = surf.knotvector_v

        # Recover info for the current patch from the lists
        controlPts = controlPts_patches[i]
        weights = [1]*len(controlPts)
        U = U_patches[i]
        sigma = sigma_patches[i]
        strain = strain_patches[i]
        if use_damage:
            damage = damage_patches[i]

        print('np.shape(sigma) = ', np.shape(sigma))


        shell['p'] = P
        shell['q'] = Q 
        shell['uKnot'] = uKnot 
        shell['vKnot'] = vKnot 
        shell['noPtsX'] = noPtsX 
        shell['noPtsY'] = noPtsY 
        # shell['weights'] = np.transpose(weights)
        shell['weights'] = weights
        shell['controlPts'] = controlPts
        shell['u'] = U 
        shell['stresses'] = sigma
        shell['strains'] = strain
        if use_damage:
            shell['damage'] = damage
        shell['nLayers'] = Nlayers
        shell['zLocs'] = zLocs


        matName = f'shell_{i}_{time}.mat'
        sio.savemat(matName,shell)

        print("Done saving .mat = ", matName)







