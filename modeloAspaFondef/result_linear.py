import scipy as sp
import numpy as np

from matplotlib.pylab import *

# 5 second of acceleration
defs_linear=loadtxt("Node3100_Y_linear85.out")
# defs_gnl=loadtxt("Node3100_Y_GNL.out")
# defs_gnl=loadtxt("Node3100_2_Y_GNL.out")

# Linear
plot(defs_linear[:,0], 100*(defs_linear[:,1]-1*defs_linear[0,1]),'b',label = 'Linear',alpha=0.7)

# # GNL
# plot(defs_gnl[:,0], 100*(defs_gnl[:,1]-1*defs_gnl[0,1]),'r',label = 'Non-Linear',alpha=0.5)

plot(defs_linear[:,0], 100*np.array([0.14/2]*len(defs_linear[:,0])),'--g')
plot(defs_linear[:,0], 100*np.array([-0.14/2]*len(defs_linear[:,0])),'--g')

legend()
ylabel('Y disp [cm]')
xlabel('Time [s]')
show()