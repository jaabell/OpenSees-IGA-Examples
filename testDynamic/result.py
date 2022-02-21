import scipy as sp

from matplotlib.pylab import *



defs=loadtxt("Node2.out")
# defs=loadtxt("Node.out")
plot(defs[:,0], defs[:,1],'-o')
show()

defs12_Z=loadtxt("Node12_Z.out")
defs12_X=loadtxt("Node12_X.out")


plot(defs12_Z[:,0], defs12_Z[:,1])
plot(defs12_X[:,0], defs12_X[:,1])
show()