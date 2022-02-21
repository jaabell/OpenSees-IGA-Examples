import scipy as sp
import numpy as np

from matplotlib.pylab import *

# 2 seconds of acceleration
defs_linear=loadtxt("Linear10HZ.out")
# defs_gnl=loadtxt("Node3100_Y_GNL.out")
defs_gnl=loadtxt("GNL10HZ.out")

# Linear
plot(defs_linear[:,0], 100*(defs_linear[:,1]-1*defs_linear[0,1]),'b',label = 'Linear',alpha=0.7)

# GNL
plot(defs_gnl[:,0], 100*(defs_gnl[:,1]-1*defs_gnl[0,1]),'r',label = 'Non-Linear',alpha=0.5)

plot(defs_linear[:,0], 100*np.array([0.14/2]*len(defs_linear[:,0])),'--g', label = 'Base displacement')
plot(defs_linear[:,0], 100*np.array([-0.14/2]*len(defs_linear[:,0])),'--g')

legend()
ylabel('Y disp [cm]')
xlabel('Time [s]')
title("10 Hz oscillation, vertical displacement at tip")

# LLEGO HASTA 0.5 ANTES DEL CAMBIO EN BENDING STRIP

# from freq import *

# amplitude = 0.14 / 2 # 14 cm de carrito
# t_steady = 3
# tMax = 4
# ω_min = 0   # Hz 
# ω_max = 10   # Hz
# nPoints_accel=700
# nPoints_steady=400
# [deltaT_vector,t,dispY] = generateOscillation(amplitude,t_steady,tMax,ω_min, ω_max,nPoints_accel,nPoints_steady)
# dispYpp=np.gradient(dispY)
# plot(t,100*dispYpp, '-r')
# plot(t,dispY,'-b')
# plot(t,100*dispY,'og')


show()