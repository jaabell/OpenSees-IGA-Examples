import scipy as sp
import numpy as np

from matplotlib.pylab import *

from scipy.interpolate import interp1d

fig=figure()

i_cycle = 0
fatigue = True


def plotAndUpdate(i):

	fig.clear()

	ax1 = fig.add_subplot(3,1,1)

	# subplot(3,1,1)

	# 2 seconds of acceleration
	# defs_linear=loadtxt("Linear8HZ_2.out")
	# defs_linear=loadtxt("Linear_fatigue.out")
	
	# defs_linear_tip=loadtxt("Linear_fatigue_tip.out")
	# defs_linear_nf=loadtxt("Linear8Hz_Notfilled.out")
	# defs_gnl_tip=loadtxt("GNL8HZ_tip.out")
	# defs_gnl=loadtxt("GNL8HZ.out")

	useTip = False

	if useTip:
		if fatigue:
			# defs_gnl = loadtxt(f"GNL_fatigue_tip_{i_cycle}.out")
			defs_linear = loadtxt(f"Linear_fatigue_tip_{i_cycle}.out")
		else:
			defs_gnl = loadtxt("GNL_fatigue_tip.out")
			defs_linear = loadtxt("Linear_fatigue_tip.out")
	else:
		if fatigue:
			# defs_gnl=loadtxt(f"GNL_fatigue_{i_cycle}.out")
			defs_linear=loadtxt(f"Linear_fatigueSimple_{i_cycle}.out")
			# defs_linear=loadtxt(f"Linear_fatigueSimple_{i_cycle}.out")
		else:
			defs_gnl=loadtxt("GNL_fatigue.out")
			defs_linear=loadtxt("Linear_fatigue.out")
	  
	# Linear


	# GNL
	# ax1.plot(defs_gnl[:,0], 100*(defs_gnl[:,1]-1*defs_gnl[0,1]),'r',label = 'GNL',alpha=0.5)


	# start_time = 1
	start_time = 0
	# input_path = '/home/rodrigo/FelipeElgueta/OpenSees-IGA-Examples/modeloAspaFondef/experimental/iga_inputs/case01_190715_13_20_05_y_tip.txt'
	input_path = '/home/rodrigo/FelipeElgueta/OpenSees-IGA-Examples/modeloAspaFondef/experimental/iga_inputs/case00_190609_11_21_47_ypp_accel.txt'

	accelY_input = np.loadtxt(input_path,skiprows=1)
	nData_input = len(accelY_input)
	dT_input = 0.00333
	t_input = np.arange(0,nData_input*dT_input,dT_input)
	accelY_input_int = interp1d(t_input,accelY_input)
	accelY_input_interpolated = accelY_input_int(defs_linear[:,0])

	g = 9.807

	t = defs_linear[:,0] + start_time
	accelY_input_interpolated = accelY_input_int(t) * g

	accel_measured_absolute = accelY_input_interpolated


	ax1.plot(t-start_time,accelY_input_interpolated, '-g', label = 'Measured acceleration')

	# # Base displacements
	# ax1.plot(defs_linear[:,0], 100*np.array([0.14/2]*len(defs_linear[:,0])),'--g',label = 'Base displacement')
	# ax1.plot(defs_linear[:,0], 100*np.array([-0.14/2]*len(defs_linear[:,0])),'--g')



	# Getting data from input file
	# start_time = 1
	start_time = 0
	# subplot(3,1,2)
	ax2 = fig.add_subplot(3,1,2)

	# input_path = '/home/rodrigo/FelipeElgueta/OpenSees-IGA-Examples/modeloAspaFondef/experimental/iga_inputs/case01_190715_13_20_05_y.txt'
	input_path = '/home/rodrigo/FelipeElgueta/OpenSees-IGA-Examples/modeloAspaFondef/experimental/iga_inputs/case00_190609_11_21_47_ypp.txt'
	accelBase_input = np.loadtxt(input_path,skiprows=1)
	nData_input = len(accelBase_input)
	dT_input = 0.00333
	t_input = np.arange(0,nData_input*dT_input,dT_input)


	accelBase_input_int = interp1d(t_input,accelBase_input)
	# dispY_input_interpolated = dispY_input_int(defs_linear[:,0])
	accelBase_input_interpolated = accelBase_input_int(t_input)

	t = defs_linear[:,0] + start_time

	# t -= start_time
	accelBase_input_interpolated = accelBase_input_int(t) * g

	accel_linear_rel = defs_linear[:,1]
	accel_linear_absolute = accel_linear_rel + accelBase_input_interpolated
	ax1.plot(defs_linear[:,0], accel_linear_absolute,'b',label = 'Linear',alpha=0.7)
	# ax1.plot(defs_linear[:,0], accel_linear_absolute,'o',alpha=0.7)


	# t = defs_gnl[:,0] + start_time
	# accelBase_input_interpolated = accelBase_input_int(t) * g
	# accel_gnl_rel = defs_gnl[:,1]
	# accel_gnl_absolute = accel_gnl_rel + 1*accelBase_input_interpolated
	# ax1.plot(defs_gnl[:,0], accel_gnl_absolute,'r',label = 'GNL',alpha=0.7)


	ax1.legend()
	# ax1.ylabel('Y disp [cm]')
	ax1.set_ylabel('Y accel [m/s2]')
	# ax1.xlabel('Time [s]')
	ax1.set_xlabel('Time [s]')
	ax1.set_title("4.43 Hz fatigue test, vertical acceleration")
	ax1.grid()



	ax2.plot(t,accelBase_input_interpolated, label = 'Base acceleration')
	# ax2.ylabel('Base disp [cm]')
	ax2.set_ylabel('Base accel [m/s2]')
	ax2.set_xlabel('Time [s]')
	ax2.legend()
	ax2.grid()

	# Base displacement history
	# subplot(3,1,3)
	ax3 = fig.add_subplot(3,1,3)

	ax3.plot(t_input,accelBase_input)
	# t_input -= start_time
	# ax3.plot(t,100*dispY_input_int(t))
	yMin = -max(abs(accelBase_input))
	ax3.vlines(defs_linear[:,0][-1]+start_time, ymin = yMin, ymax = -yMin, linestyles = 'dashed', colors = 'red', label = 'elapsed time linear')
	# ax3.vlines(defs_gnl[:,0][-1]+start_time, ymin = yMin, ymax = -yMin, linestyles = 'dashed', colors = 'green', label = 'elapsed time GNL')

	ax3.set_ylabel('Base accel \n history [m/s2]')
	ax3.set_xlabel('Time [s]')
	ax3.legend()

	tight_layout()
	ax3.grid()


	# show()


import matplotlib.animation as animation

ani = animation.FuncAnimation(fig,plotAndUpdate,interval=10*1000)
show()