import scipy as sp
import numpy as np

from matplotlib.pylab import *

# 2 seconds of acceleration
defs_linear=loadtxt("Linear16HZ.out")
# defs_gnl=loadtxt("Node3100_Y_GNL.out")
# defs_gnl=loadtxt("GNL8HZ.out")
  
# Linear
plot(defs_linear[:,0], 100*(defs_linear[:,1]-1*defs_linear[0,1]),'b',label = 'Linear',alpha=0.7)

# GNL
# plot(defs_gnl[:,0], 100*(defs_gnl[:,1]-1*defs_gnl[0,1]),'r',label = 'Non-Linear',alpha=0.5)

plot(defs_linear[:,0], 100*np.array([0.14/2]*len(defs_linear[:,0])),'--g',label = 'Base displacement')
plot(defs_linear[:,0], 100*np.array([-0.14/2]*len(defs_linear[:,0])),'--g')

legend()
ylabel('Y disp [cm]')
xlabel('Time [s]')
title("16 Hz oscillation, vertical displacement at tip")


# LLEGO HASTA 1.83 ANTES con newmark !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# from freq import *

# amplitude = 0.14 / 2 # 14 cm de carrito
# t_steady = 2
# tMax = 4
# ω_min = 0   # Hz 
# ω_max = 8   # Hz
# nPoints_accel=300		
# nPoints_steady=300
# [deltaT_vector,t,dispY] = generateOscillation(amplitude,t_steady,tMax,ω_min, ω_max,nPoints_accel,nPoints_steady)
# dispYpp=np.gradient(dispY)
# plot(t,100*dispYpp, '-r')
# plot(t,dispY,'-b')
# plot(t,100*dispY,'og')



# show()
# print(deltaT_vector)
# print(t)

show()