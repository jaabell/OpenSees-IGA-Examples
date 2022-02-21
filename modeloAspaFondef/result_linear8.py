import scipy as sp
import numpy as np

from matplotlib.pylab import *



defs_linear=loadtxt("Node3100_Y_linear9.out")


# Linear
plot(defs_linear[:,0], 100*(defs_linear[:,1]-1*defs_linear[0,1]),'-b',label = 'Linear',alpha=0.7)


plot(defs_linear[:,0], 100*np.array([0.14/2]*len(defs_linear[:,0])),'--g')
plot(defs_linear[:,0], 100*np.array([-0.14/2]*len(defs_linear[:,0])),'--g')

legend()
ylabel('Y disp [cm]')
xlabel('Time [s]')
show()



