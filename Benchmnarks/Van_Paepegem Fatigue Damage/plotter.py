import numpy as np
import scipy as sp

data=sp.loadtxt("Node12_Z.out")

loadFactor=data[:,0]
uTip=data[:,1]

from matplotlib.pylab import *

plot(uTip,loadFactor)
show()
plot(loadFactor)
show()
