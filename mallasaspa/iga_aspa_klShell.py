# Cantilever beam verification using Isogeometric Kirchhoff Love Shell
#  python using geomdl
# Felipe Elgueta
# Universidad de los Andes, Chile

import os
import sys
# from geomdl import BSpline
from geomdl import NURBS
from geomdl import operations
from geomdl import exchange
from geomdl import utilities
from geomdl import helpers
from geomdl import compatibility
from geomdl.visualization import VisVTK as vis

from scipy.sparse import csr_matrix, lil_matrix, linalg
from scipy.sparse.linalg import eigsh

from fileGetter import surfFromFile, generateKnotVector

surf=surfFromFile("coarse_1.dat",deg_elevate=True)



# Fix file path
os.chdir(os.path.dirname(os.path.realpath(__file__)))

# add path to other files
sys.path.append("../data/")
sys.path.append("../visualization/")
sys.path.append("../nurbs-util/")
sys.path.append("../fem-functions/")
sys.path.append("../meshing/")
sys.path.append("../fem_util/")
sys.path.append("../nurbs-util/")
sys.path.append("../python_util/")
sys.path.append("../examples/")
sys.path.append("../C_files_win/argout/")


from surfVisualize import surfVisualize
from getCtrlPtsAndWeights import *
from elasticityMatrix import elasticityMatrix
from generateIGA2DMesh import generateIGA2DMesh
from quadrature import quadrature
from parent2ParametricSpace import *
from jacobianPaPaMapping import *
from fill_input import *
from argout import *


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

# Create a BSpline surface instance
# surf = NURBS.Surface()



# import geometry data
# from iga_cantilever_klShell_data_test import *

# Set degrees
# surf.degree_u = p-1
# surf.degree_v = q-1

# Set control points
# surf.set_ctrlpts(controlPts, 4, 4)

# Set knot vectors

# surf.knotvector_u=generateKnotVector(surf.degree_u, 4)
# surf.knotvector_v=generateKnotVector(surf.degree_u, 4)

# Refine knot vectors and update geometry (refine both directions)
                                #u,v
refU=0
refV=0
operations.refine_knotvector(surf, [refU, refV])


# Get new control points and weights after refining
weights = surf.weights
controlPts = surf.ctrlpts2d[:]

# Update geometry after refining
noPtsX = surf.ctrlpts_size_u
noPtsY = surf.ctrlpts_size_v

# Flip control points to [u][v] (they are given by the getters in [v][u])
controlPts, weights = getCtrlPtsAndWeights(surf)
# controlPts = surf.ctrlpts


# # Set degrees
# if refU==0:
#     p=2
# if refV==0:
#     q=2
# surf.degree_u = p
# surf.degree_v = q

# uKnot=generateKnotVector(surf.degree_u, surf.ctrlpts_size_u)
# surf.knotvector_u=uKnot
# vKnot=generateKnotVector(surf.degree_v, surf.ctrlpts_size_v)
# surf.knotvector_v=vKnot

uKnot=surf.knotvector_u
vKnot=surf.knotvector_v

p=surf.degree_u
q=surf.degree_v

noPtsX = surf.ctrlpts_size_u
noPtsY = surf.ctrlpts_size_v

noGPs = surf.degree_u + 1

noCtrPts = len(controlPts)
noDofs = noCtrPts * 3 # three displacement dofs per node

#Visualize surface
surfVisualize(surf, hold=True)
# exit()


# Material parameters
t=1
E = 2.03e11  # Young's modulus N/m^2
nu = 0.3  # Poisson's ratio
rho = 7.7e03 #*9.807 # kg/m^3 

memStiff = E * t / (1 - nu**2)
benStiff = E * t**3 / 12 / (1 - nu**2)

# find boundary nodes for boundary conditions
leftNodes=(np.arange(0, noPtsX))*(noPtsY)
leftNodes=np.arange(noPtsX)
fixedNodes=leftNodes

# print("leftNodes: ", leftNodes)
# # exit()

# nextToLeftNodes=np.arange(noPtsX, 2*noPtsX)
# fixedNodes=np.append(fixedNodes,nextToLeftNodes)

print("fixedNodes: ", fixedNodes)
# exit()
#Find free nodes
freeNodes=list(set(np.arange(noCtrPts))-set(fixedNodes))

# % essential boundary conditions (fixed along edge)
udofs      = 3*(fixedNodes+1)-2-1
vdofs      = 3*(fixedNodes+1)-1-1
wdofs      = 3*(fixedNodes+1)-0-1



uFixed = np.zeros((len(fixedNodes)))
vFixed = np.zeros((len(fixedNodes)))
wFixed = np.zeros((len(fixedNodes)))

oldDofs=range(noDofs)
freeDofs = (list(set(oldDofs) - set(udofs)-set(vdofs)-set(wdofs)))


# build connectivity ...
noElems, noElemsU, noElemsV, index, elRangeU, elRangeV, element, elConnU, elConnV = generateIGA2DMesh(
    surf)

# Gauss quadrature rule
noGPs = 6
noGpEle = noGPs**2
[W, Q] = quadrature(noGPs, 'GAUSS', 2)  # noGPs x noGPs point quadrature


# % initialization

K = lil_matrix((noDofs, noDofs))
M = lil_matrix((noDofs, noDofs))
u = np.zeros((noDofs, 1))

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%
# %%% PROCESSING
# %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

print("ASSEMBLING THE SYSTEM")
# % Loop over elements (knot spans)

for e in range(noElems):
    idu = index[e, 0]
    idv = index[e, 1]

    xiE = elRangeU[idu, :]  # % [xi_i,xi_i+1]
    etaE = elRangeV[idv, :]  # % [eta_j,eta_j+1]

    sctr = np.array(element[e, :], dtype=int)
    nn = len(sctr)

    nn3 = 3 * nn
    sctrB = np.zeros((nn3))


    sctrB[np.ix_(np.arange(0, nn3 / 3) * 3 + 0)] = 3 * (sctr + 1) - 2 - 1
    sctrB[np.ix_(np.arange(0, nn3 / 3) * 3 + 1)] = 3 * (sctr + 1) - 1 - 1
    sctrB[np.ix_(np.arange(0, nn3 / 3) * 3 + 2)] = 3 * (sctr + 1) - 1

    pts = np.matrix((np.array(controlPts)[np.ix_(sctr)])[:, 0:3])

    for gp in range(len(W)):
        pt = Q[gp, :]
        wt = W[gp][0]

        # % compute coords in parameter space
        xi = parent2ParametricSpace(xiE, pt[0])
        eta = parent2ParametricSpace(etaE, pt[1])
        J2 = jacobianPaPaMapping(xiE, etaE)

        # fill input arguments for c function wrapped with swig
        Xi = [xi, eta]
        arguments = [Xi, uKnot, vKnot, weights]
        input_c = fill_input(arguments)

        # Get derivatives and slice the array
        R, dRdxi = Nurbs2DBasis2ndDers1(input_c, (p, q))
        dRdeta, dR2dxi = Nurbs2DBasis2ndDers2(input_c, (p, q))
        dR2det, dR2dxe = Nurbs2DBasis2ndDers3(input_c, (p, q))


        noFuncs = (p + 1) * (q + 1)
        R = R[0, 0:noFuncs]
        dRdxi = dRdxi[0, 0:noFuncs]
        dRdeta = dRdeta[0, 0:noFuncs]
        dR2dxi = dR2dxi[0, 0:noFuncs]
        dR2det = dR2det[0, 0:noFuncs]
        dR2dxe = dR2dxe[0, 0:noFuncs]

 

        #  compute the jacobian of physical and parameter domain mapping
        # then the derivative w.r.t spatial physical coordinates

        jacob = np.matrix([dRdxi, dRdeta]) * pts  # 2x2 matrix
        jacob2 = np.matrix([dR2dxi, dR2det, dR2dxe]) * pts  # 3x2 matrix

        dxdxi = jacob[0, 0]
        dydxi = jacob[0, 1]
        dxdet = jacob[1, 0]
        dydet = jacob[1, 1]

        j33 = np.matrix([[dxdxi**2, dydxi**2, 2 * dxdxi * dydxi],
                         [dxdet**2, dydet**2, 2 * dxdet * dydet],
                         [dxdxi * dxdet, dydxi * dydet, dxdxi * dydet + dxdet * dydxi]])

        # a1, a2 and a3 vectors (surface basis vectors)
        # and its derivatives

        a1 = jacob[0, :]
        a2 = jacob[1, :]
        a3 = np.cross(a1, a2)
        norma = np.linalg.norm(a3)
        a3 = a3 / norma
        J1 = norma

        a11 = jacob2[0, :]
        a22 = jacob2[1, :]
        a12 = jacob2[2, :]

        # dot products of ai and ei

        a1e1 = a1.T[0]
        a1e2 = a1.T[1]
        a1e3 = a1.T[2]
        a2e1 = a2.T[0]
        a2e2 = a2.T[1]
        a2e3 = a2.T[2]

        # R_I,2*a1 + R_I,1*a2 for all shape functions
        noBasis = len(R)
        dRIa = np.zeros((3, noBasis))
        for i in range(noBasis):
            dRIa[:, i] = dRdeta[i] * a1 + dRdxi[i] * a2

        # compute the constitutive matrix C
        a_11 = np.dot(a1, a1.T)[0, 0]  # to get the scalar
        a_12 = np.dot(a1, a2.T)[0, 0]
        a_21 = np.dot(a2, a1.T)[0, 0]
        a_22 = np.dot(a2, a2.T)[0, 0]

        aa1 = np.linalg.solve(np.matrix([[a_11, a_21],
                                         [a_12, a_22]
                                         ]), [1, 0])
        aa2 = np.linalg.solve(np.matrix([[a_11, a_21],
                                         [a_12, a_22]
                                         ]), [0, 1])

        au11 = aa1[0]
        au12 = aa1[1]
        au22 = aa2[1]

        C = np.matrix([[au11**2, nu * au11 * au22 + (1 - nu) * au12**2, au11 * au12],
                       [nu * au11 * au22 + (1 - nu) * au12 **2, au22**2, au22 * au12],
                       [au11 * au12, au22 * au12, 0.5 *((1 - nu) * au11 * au22 + (1 + nu) * au12**2)]
                       ])


        # membrane and bending B matrices
        Bmem = np.zeros((3, noBasis * 3))
        Bben = np.zeros((3, noBasis * 3))

        for i in range(noBasis):
            dRIdx = dRdxi[i]
            dRIdy = dRdeta[i]

            id_ = np.arange((i - 1 + 1) * 3 + 1 - 1, 3 * (i + 1))


            Bmem[:, id_] = [
                [dRIdx * a1e1[0, 0], dRIdx * a1e2[0, 0], dRIdx * a1e3[0, 0]],
                [dRIdy * a2e1[0, 0], dRIdy * a2e2[0, 0], dRIdy * a2e3[0, 0]],
                [dRIa[0, i],  dRIa[1, i],  dRIa[2, i]]
            ]

            BI1 = -dR2dxi[i] * a3 + 1 / norma * (dRIdx * np.cross(a11, a2) + dRIdy * np.cross(
                a1, a11) + np.dot(a3, a11.T) * (dRIdx * np.cross(a2, a3) + dRIdy * np.cross(a3, a3)))
            BI2 = -dR2det[i] * a3 + 1 / norma * (dRIdx * np.cross(a22, a2) + dRIdy * np.cross(
                a1, a22) + np.dot(a3, a22.T) * (dRIdx * np.cross(a2, a3) + dRIdy * np.cross(a3, a3)))
            BI3 = -dR2dxe[i] * a3 + 1 / norma * (dRIdx * np.cross(a12, a2) + dRIdy * np.cross(
                a1, a12) + np.dot(a3, a12.T) * (dRIdx * np.cross(a2, a3) + dRIdy * np.cross(a3, a3)))


            Bben[:, id_] = [BI1.tolist()[0], BI2.tolist()[0],
                            (2 * BI3).tolist()[0]]

            

        # compute elementary stiffness matrix and
        # assemble it to the global matrix

        # Forming N matrix
        N = np.zeros((3, 3 * noFuncs))
        for i in range(noFuncs):
            N[0, 3 * i] = R[i]
            N[1, 3 * i+1] = R[i]
            N[2, 3 * i+2] = R[i]
        # End forming N matrix

        Ke = memStiff * np.matmul(np.matmul(Bmem.T, C), Bmem) * J1 * J2 * \
            wt + benStiff * \
            np.matmul(np.matmul(Bben.T, C), Bben) * J1 * J2 * wt


        K[np.ix_(sctrB, sctrB)] = K[np.ix_(sctrB, sctrB)] + Ke

        M[np.ix_(sctrB, sctrB)] = M[np.ix_(sctrB, sctrB)] + \
            np.matmul(N.T, N) * rho * J1 * J2 * wt*t


print("APPLYING BOUNDARY CONDITIONS")

import scipy as sp
# K=K.toarray()
# bcwt = np.mean(K.diagonal())
# K[np.ix_(udofs),:]=0 #% zero out the rows and  columns of the K matrix
# K[np.ix_(vdofs),:]=0
# K[np.ix_(wdofs),:]=0
# K[:,np.ix_(udofs)]=0 
# K[:,np.ix_(vdofs)]=0
# K[:,np.ix_(wdofs)]=0
# K[np.ix_(udofs,udofs)]=bcwt*sp.eye(len(udofs))
# K[np.ix_(vdofs,vdofs)]=bcwt*sp.eye(len(vdofs))
# K[np.ix_(wdofs,wdofs)]=bcwt*sp.eye(len(wdofs))
# K=sp.sparse.csr_matrix(K)




# M=M.toarray()
# bcwt = np.mean(M.diagonal())
# M[np.ix_(udofs),:]=0 #% zero out the rows and  columns of the K matrix
# M[np.ix_(vdofs),:]=0
# M[np.ix_(wdofs),:]=0
# M[:,np.ix_(udofs)]=0 
# M[:,np.ix_(vdofs)]=0
# M[:,np.ix_(wdofs)]=0
# M[np.ix_(udofs,udofs)]=bcwt*sp.eye(len(udofs))
# M[np.ix_(vdofs,vdofs)]=bcwt*sp.eye(len(vdofs))
# M[np.ix_(wdofs,wdofs)]=bcwt*sp.eye(len(wdofs))
# M=sp.sparse.csr_matrix(M)



Kff=K[np.ix_(freeDofs,freeDofs)].tocsr()
Mff=M[np.ix_(freeDofs,freeDofs)].tocsr()

from matplotlib.pylab import *
# figure(1)
# imshow(K.todense(), cmap=cm.coolwarm)
# colorbar()
# show()

# figure(1)
# imshow(M.todense(), cmap=cm.coolwarm)
# colorbar()
# show()


def modal(nModes,save=False):

    print("GETTING MODAL SHAPES")

    w, phi_full = eigsh(Kff, M=Mff,k=60,sigma=0, maxiter=1000000000)

    from scipy import linalg

    w=sp.sqrt(w)
    order=argsort(w)
    w=w[order]
    T=2*sp.pi/w.real
    # for (item,i) in zip(T,range(10)):
    #     print('{:f}'.format(item))

    T_real=[0.0179,0.00732,0.00292,0.00228,0.00201]
    print("Comparing periods with reference")
    for i in range(5):
        print("T[i]: ", T[i])
        print("T_real[i]: ", T_real[i])
        print("T[i]/T_real[i]: ", T[i]/T_real[i])
        print("\n")

    for mode in range(nModes):
        nMode=mode
        phi=phi_full[:,order[nMode]]

        #Getting control pts from surface
        controlPts_1=surf.ctrlpts2d[:]
        controlPts_1=compatibility.flip_ctrlpts2d(surf.ctrlpts2d[:])
        controlPts_1 = (np.array(controlPts_1).reshape(surf.ctrlpts_size_u * surf.ctrlpts_size_v, 4)).tolist()

        #Add modal shape deformation to control points
        fac=10*t
        for (i,k) in zip(freeNodes,range(3*len(freeNodes))):
            controlPts_1[i][0]+=fac*phi[3*k]
            controlPts_1[i][1]+=fac*phi[3*k+1]
            controlPts_1[i][2]+=fac*phi[3*k+2]


        # Create a new deformes Nurbs surface
        surf_1 = NURBS.Surface()

        # Set degrees
        surf_1.degree_u = p
        surf_1.degree_v = q

        # Set control points
        surf_1.set_ctrlpts(controlPts_1, surf.ctrlpts_size_v, surf.ctrlpts_size_u)





        # Set knot vectors
        surf_1.knotvector_u=generateKnotVector(surf.degree_u, surf_1.ctrlpts_size_u)
        surf_1.knotvector_v=generateKnotVector(surf.degree_v, surf_1.ctrlpts_size_v)

        # surfVisualize(surf, hold=True)
        surfVisualize(surf_1,mode,save=save)

modal(10)





