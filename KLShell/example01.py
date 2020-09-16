import numpy as np
import opensees as ops

La = 10.0  	#
Lb = 10.0  	#
t = 0.05  	# m


ops.model('basic', '-ndm', 3, '-ndf', 3)


uKnot = np.array([ 0.0,0.0,0.0,0.5,1.0,1.0,1.0,])
vKnot = np.array([ 0.0,0.0,0.0,0.5,1.0,1.0,1.0,])
controlPts = np.array([
    [-La / 2, -Lb / 2, 0, 1], 
    [-La / 2, -Lb / 2 / 2, 0,1], 
    [-La / 2, Lb / 2 / 2, 0, 1], 
    [-La / 2, Lb / 2, 0, 1],
    [-La / 2 / 2, -Lb / 2, 0, 1], 
    [-La / 2 / 2, -Lb / 2 / 2, 0,1], 
    [-La / 2 / 2, Lb / 2 / 2, 0, 1], 
    [-La / 2 / 2, Lb / 2, 0, 1],
    [La / 2 / 2, -Lb / 2, 0, 1], 
    [La / 2 / 2, -Lb / 2 / 2, 0,1], 
    [La / 2 / 2, Lb / 2 / 2, 0, 1], 
    [La / 2 / 2, Lb / 2, 0, 1],
    [La / 2, -Lb / 2, 0, 1], 
    [La / 2, -Lb / 2 / 2, 0,1], 
    [La / 2, Lb / 2 / 2, 0, 1], 
    [La / 2, Lb / 2, 0, 1]
])

patchTag = 1
P = 3
Q = 3


ops.IGA("Patch",patchTag, P, Q,"-uKnot", *uKnot, "-vKnot", *vKnot, "-controlPts", *controlPts.flatten())
