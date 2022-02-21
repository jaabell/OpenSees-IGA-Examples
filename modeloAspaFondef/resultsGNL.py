import scipy as sp
import numpy as np

from matplotlib.pylab import *

# 5 seconds of acceleration
# defs_gnl=loadtxt("Node3100_Y_GNL.out")
defs_gnl=loadtxt("Node3100_2_Y_GNL.out")


# GNL
plot(defs_gnl[:,0], 100*(defs_gnl[:,1]-1*defs_gnl[0,1]),'r',label = 'Non-Linear',alpha=0.5)

plot(defs_gnl[:,0], 100*np.array([0.14/2]*len(defs_gnl[:,0])),'--g')
plot(defs_gnl[:,0], 100*np.array([-0.14/2]*len(defs_gnl[:,0])),'--g')

legend()
ylabel('Y disp [cm]')
xlabel('Time [s]')
show()