from scipy.sparse.linalg import eigsh
from geomdl.visualization import VisVTK, VisMPL
from geomdl import exchange_vtk
from geomdl import exchange
from edgeHandler import *
from scipy import interpolate
import numpy as np
import opensees as ops
from math import *




from six.moves import cPickle as pickle #for performance

def load_dict(filename_):
    with open(filename_, 'rb') as f:
        ret_di = pickle.load(f)
    return ret_di


surfaces_data_loaded = load_dict('surfaces_data.pkl')

surfList=[]

for surface_name in surfaces_data_loaded:
    surface=surfaces_data_loaded[surface_name]
    # for key in surface:
    #     print (key, surface[key])
    # print('\n')



    nPtsU = surface['nPtsU']
    nPtsV = surface['nPtsV']
    uKnot = surface['uKnot']
    vKnot = surface['vKnot']
    p = surface["p"]
    q = surface["q"]


    nodeStartTag= surface['nodeStartTag']
    nodeEndTag= surface['nodeEndTag']

    nodeTags=np.arange(nodeStartTag,nodeEndTag+1)
    nodeTags=nodeTags.reshape([nPtsU,nPtsV],order = 'F')

    surface["nodeTags"] = nodeTags

    # print('nodeTags = \n', nodeTags, '\n')

    # controlPts = surface['ctrlpts'].tolist()
    controlPts = surface['ctrlpts'].reshape([nPtsU * nPtsV, 4]).tolist()

    # surf.set_ctrlpts(controlPts.reshape(nPoints, 4).tolist(), surf.ctrlpts_size_u, surf.ctrlpts_size_v)


    # Create a NURBS surface instance
    surf = NURBS.Surface()

    # Set degrees
    surf.degree_u = p
    surf.degree_v = q

    # Set control points
    surf.set_ctrlpts(controlPts, nPtsU, nPtsV)

    nPoints = surf.ctrlpts_size_u * surf.ctrlpts_size_v
    shape = np.array(compatibility.flip_ctrlpts2d(surf.ctrlpts2d[:])).shape
    controlPts = np.array(controlPts).reshape(shape)
    controlPts = np.array(compatibility.flip_ctrlpts2d(controlPts))
    surf.set_ctrlpts(controlPts.reshape(nPoints, 4).tolist(), surf.ctrlpts_size_u, surf.ctrlpts_size_v)

    # Set knot vectors
    # surf.knotvector_u = knotvector.generate(surf.degree_u, surf.ctrlpts_size_u)
    # surf.knotvector_v = knotvector.generate(surf.degree_v, surf.ctrlpts_size_v)

    surf.knotvector_u = uKnot
    surf.knotvector_v = vKnot

    surfType = surface['type']

    if surfType == 'BendingStrip':
        continue


    surfList.append(surf)

container = multi.SurfaceContainer(surfList)
container.sample_size = 20

for surf in container:
    surf.evaluate()


# Visualization configuration
container.vis = VisVTK.VisSurface(ctrlpts=False)

# Render the aspa
container.render()



# PostProcessing

import scipy.io as sio
import h5py

import re

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]

use_damage = True


import glob 

# Full model
# fileNames = sorted(list(glob.glob("/home/rodrigo/FelipeElgueta/OpenSees-IGA-Examples/modeloAspaFondef/4_43Hz_accel/aspa_linear_*.mpco")))
# shellName = 'shell'

# Simple model
fileNames = sorted(list(glob.glob("/home/rodrigo/FelipeElgueta/OpenSees-IGA-Examples/modeloAspaFondef/calibrationCases/case_0/aspa_linearSimple_*.mpco")))  # Simple model, just to check the damage distribution is continous
shellName = 'shellSimple'

# fileNames.sort(key=lambda f: int(''.join(filter(str.isdigit, f))))
fileNames.sort(key=natural_keys)

print('len(fileNames) = ', len(fileNames))
print("fileNames = ", fileNames)

for file in fileNames:

    file = h5py.File(file, 'r')
    # stages = [4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34]

    # Getting stages info
    fileInfo = file
    stages = [name for name in fileInfo.keys()][1:] # [1:] is to remove the first stage that is info
    # stages.sort()
    stages.sort(key=natural_keys)

    print("stages = ", stages)
    print("len(stages) = ", len(stages))


    startTime = 1
    
    maxTimes = 1

    currTime = -1

    for stage in stages:

        currTime += 1

        modelStage = stage

        patches = file[modelStage]['MODEL']['ELEMENTS']['1-IGASurfacePatch[0:0]'][:, 0] # These are the subdomains -> Iga Patches
        elements = file[modelStage]['MODEL']['ELEMENTS']['253-IGAKLShell[0:0]'][:, 0] # These are the IGA_KLShell elements

        nodes = file[modelStage]['MODEL']['ELEMENTS']['253-IGAKLShell[0:0]'][:, 1:]

        numPatches = 1

        nodesOnPatches = []
        elementsOnPatches = []

        continueSearch = True
        while continueSearch:

          for i in range(1,len(nodes)):
            currentElement = elements[i-1]
            nextElement = elements[i]

            # print("currentElement = ", currentElement)
            # print("nextElement = ", nextElement)

            foundPatch = True

            if abs(nextElement - currentElement) > 1:
              continueSearch = True
              foundPatch = True
              numPatches += 1
              # print('breaking in i = ', i)
              break_i = i
              nodesOnPatches.append(nodes[:i,:])
              elementsOnPatches.append(elements[:i])
              nodes = nodes[i:,:]
              elements = elements[i:]
              # print('nodes = ', nodes)
              # print('elements = ', elements)
              break
            else:
              foundPatch = False

          if not foundPatch:
            continueSearch = False
            nodesOnPatches.append(nodes)
            elementsOnPatches.append(elements)

          else:
            pass

        # Getting layup info
        layupData = file[modelStage]['MODEL']['SECTION_ASSIGNMENTS']
        layupData_keys = layupData.keys()
        numLayups = len(layupData)
        # print('numLayups = ', numLayups)
        # print('layupData_keys = ', layupData_keys)

        zLocs_list = []

        for key in layupData_keys:
            # print('key = ', key)

            layupData_here = layupData[key]

            zLocs_layups = layupData_here['FIBER_DATA'][:,1] # Shape (numMaterial x 3 (algo, zLoc, thickness))

            # print('zLocs_layups = ', zLocs_layups)

            zLocs_list.append(zLocs_layups)




        # Getting nodeLimits in order to know which nodes are from which patch
        nodeLimits = []
        for nodes in nodesOnPatches:
          nodeLimits.append([nodes[0,0],nodes[-1,-1]])

        nTimes = len(file[modelStage]['RESULTS']['ON_NODES']['DISPLACEMENT']['DATA'])

        print('nTimes = ', nTimes)

        print('stage = ', stage)

        firstTime = int(list(file[modelStage]['RESULTS']['ON_NODES']['DISPLACEMENT']['DATA'].keys())[0].strip('STEP_'))
        print(firstTime)
        startTime = firstTime
        print('startTime = ', startTime)

        for time in range(nTimes):

          step = f'STEP_{time + startTime}' 

          # print('step = ', step)

          controlPts = file[modelStage]['MODEL']['NODES']['COORDINATES'][:, :] # All the control points / nodes in the model
          U = file[modelStage]['RESULTS']['ON_NODES']['DISPLACEMENT']['DATA'][step][:,:] # All the displacements in the model

          nResults = len(file[modelStage]['RESULTS']['ON_ELEMENTS']['section.fiber.stress']) # Doesn't have to be stress, just need the number of files

          # The files are different because they are not the same size (nLayers)
         
          dictResults = {}

          for nFile in range(nResults):
            # Getting current file tag

            actualFile = f'253-IGAKLShell[0:0:{nFile}]'


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
            # print('actualFile = ', actualFile)


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

            # print('elementsOnPatches = ', elementsOnPatches)

            # Have to do a new grid
            if newFile:
              # I'm in the first patch
              # print("New file!")
              numEleThisPatch = len(elementsThisPatch)
              id_elements = np.arange(0,numEleThisPatch,1)

            else:
              # print('Same file!')
              elementsPrevPatch = elementsOnPatches[i-1]
              # print("elementsPrevPatch = ", elementsPrevPatch)
              # print("elementsThisPatch = ", elementsThisPatch)
              numEleThisPatch = len(elementsThisPatch)
              numElePrevPatch = len(elementsPrevPatch)
              # print('numElePrevPatch = ', numElePrevPatch)
              id_elements = np.arange(id_elements[-1]+1, id_elements[-1]+numEleThisPatch+1,1)

            # print('id_elements = ', id_elements)
            # print('len(id_elements) = ', len(id_elements))


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

            # print('id_elements = ',id_elements)
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


            # print('sigmaXX = ', sigmaXX)

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
              # print('np.shape(sigmaXX) = ',np.shape(sigmaXX))

              # Getting stresses and strains corresponding to actual layer 
              sigmaXX_layer = sigmaXX[:,:,np.ix_(rangeLayer)]
              sigmaYY_layer = sigmaYY[:,:,np.ix_(rangeLayer)]
              sigmaXY_layer = sigmaXY[:,:,np.ix_(rangeLayer)]

              # print("sigmaXX_layer = ", sigmaXX_layer)


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

            print('i = ', i)
            print("len(surfList) = ", len(surfList))
            print("Nlayers_patches = ", Nlayers_patches)
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
            shell['done'] = False

            
            print("Time = ", time)

            matName = f'mats/{shellName}_{i}_{time+currTime}.mat'

            # print('np.shape(controlPts) = ', np.shape(controlPts))
            # print('noPtsX = ', noPtsX)
            # print('noPtsY = ', noPtsY)


            try:
                print("Trying")
                shell = sio.loadmat(matName)
                print("shell['done'] = ", shell['done'])
                if shell['done']:
                    print("shell['done'] = ", shell['done'])
                    print(".mat Already done = ", matName)
                    continue
                else:
                    print("Doing")
                    pass
            except:
                pass

            sio.savemat(matName,shell)
            print("Done saving .mat = ", matName)

            # exit()

            # print('Nlayers = ', Nlayers)
            # print('len(zLocs) = ', len(zLocs))