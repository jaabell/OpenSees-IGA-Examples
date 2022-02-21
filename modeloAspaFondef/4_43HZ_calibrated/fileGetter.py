import os
from geomdl import NURBS, compatibility, knotvector, operations
import numpy as np

# Fix file path
os.chdir(os.path.dirname(os.path.realpath(__file__)))

def generateKnotVector(deg,nPts):

    import numpy as np
    knotVector=np.zeros(nPts+deg+1)
    nMiddle=len(knotVector)-2*(deg+1)
    step=1.0/(nMiddle+1)

    ini=np.zeros(deg+1)
    if step==0.5:
        middle=np.array([step])
    else:
        middle=np.arange(0+step, 1-step, step)
        middle=np.linspace(0+step, 1-step,nPts+deg+1-2*(deg+1))
    fin=ini+1

    knotVector=np.copy(ini)
    knotVector=np.append(knotVector,middle)
    knotVector=np.append(knotVector,fin)


    return knotVector


def surfFromFile(fileName, refine_U = False):

    file = open(fileName, 'r')

    nDeg = int(file.readline().strip())

    deg_u, deg_v = [int(el) for el in file.readline().strip().split()]

    noPtsV, noPtsU = [int(el) for el in file.readline().strip().split()]

    knotU = [float(el) for el in file.readline().strip().split()]
    knotV = [float(el) for el in file.readline().strip().split()]

    controlPts = []
    line = file.readline()
    while line != '':
        pt = [float(el) for el in line.strip().split()]
        controlPts.append(pt)
        line = file.readline()
    controlPts.pop(-1)
    file.close()

    # Create a NURBS surface instance
    surf = NURBS.Surface()

    # Set degrees
    surf.degree_u = deg_u 
    surf.degree_v = deg_v 

    # Set control points
    surf.set_ctrlpts(controlPts, noPtsU, noPtsV)

    # Set knot vectors
    surf.knotvector_u = knotvector.generate(surf.degree_u, surf.ctrlpts_size_u)
    surf.knotvector_v = knotvector.generate(surf.degree_v, surf.ctrlpts_size_v)

    # surf.knotvector_u = knotV
    # surf.knotvector_v = knotU


    # Refino
    if refine_U:
        operations.refine_knotvector(surf,[1,0])
    
    # test=compatibility.flip_ctrlpts2d(surf.ctrlpts2d)
    # # test=surf.ctrlpts2d[:]
    # test=np.array(test).reshape(noPtsU*noPtsV,4)
    # test=test.tolist()

    # surf.set_ctrlpts(test, noPtsU, noPtsV)



    # surf.ctrlpts_size_u = noPtsU
    # surf.ctrlpts_size_v = noPtsV



    # surf.ctrlpts=controlPts


    return surf
