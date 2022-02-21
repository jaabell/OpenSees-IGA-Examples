import scipy as sp
import numpy as np

from matplotlib.pylab import *

# 2 seconds of acceleration
# defs_linear=loadtxt("Linear8HZ_2.out")
defs_linear=loadtxt("Linear8Hz_filled.out")
defs_linear_tip=loadtxt("Linear8Hz_filled_tip.out")
# defs_linear_nf=loadtxt("Linear8Hz_Notfilled.out")
defs_gnl_tip=loadtxt("GNL8HZ_tip.out")
defs_gnl=loadtxt("GNL8HZ.out")

useTip = True

if useTip:
	defs_gnl = loadtxt("GNL8HZ_tip.out")
	defs_linear = loadtxt("Linear8Hz_filled_tip.out")
else:
	defs_gnl=loadtxt("GNL8HZ.out")
	defs_linear=loadtxt("Linear8Hz_filled.out")
  
# Linear
plot(defs_linear[:,0], 100*(defs_linear[:,1]-1*defs_linear[0,1]),'b',label = 'Linear',alpha=0.7)

# GNL
plot(defs_gnl[:,0], 100*(defs_gnl[:,1]-1*defs_gnl[0,1]),'r',label = 'Geometrically Non-Linear',alpha=0.5)

# Base displacements
plot(defs_linear[:,0], 100*np.array([0.14/2]*len(defs_linear[:,0])),'--g',label = 'Base displacement')
plot(defs_linear[:,0], 100*np.array([-0.14/2]*len(defs_linear[:,0])),'--g')

legend()
ylabel('Y disp [cm]')
xlabel('Time [s]')
title("8 Hz oscillation, vertical displacement at tip")


# LLEGO HASTA 1.9872973887034735 ANTES con newmark y con damping y factor 2e1 y espesor 1t !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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