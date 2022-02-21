import scipy as sp
import numpy as np

from matplotlib.pylab import *

from scipy.interpolate import interp1d

fig=figure()

i_cycle = 2
fatigue = True


def plotAndUpdate(i):

	fig.clear()

	ax1 = fig.add_subplot(3,1,1)

	# subplot(3,1,1)

	# 2 seconds of acceleration
	# defs_linear=loadtxt("Linear8HZ_2.out")
	defs_linear=loadtxt("Linear_fatigue.out")
	defs_linear_tip=loadtxt("Linear_fatigue_tip.out")
	# defs_linear_nf=loadtxt("Linear8Hz_Notfilled.out")
	# defs_gnl_tip=loadtxt("GNL8HZ_tip.out")
	# defs_gnl=loadtxt("GNL8HZ.out")

	useTip = False

	if useTip:
		if fatigue:
			defs_gnl = loadtxt(f"GNL_fatigue_tip_{i_cycle}.out")
			defs_linear = loadtxt(f"Linear_fatigue_tip_{i_cycle}.out")
		else:
			defs_gnl = loadtxt("GNL_fatigue_tip.out")
			defs_linear = loadtxt("Linear_fatigue_tip.out")
	else:
		if fatigue:
			defs_gnl=loadtxt(f"GNL_fatigue_{i_cycle}.out")
			# defs_linear=loadtxt(f"Linear_fatigue_{i_cycle}.out")
		else:
			defs_gnl=loadtxt("GNL_fatigue.out")
			defs_linear=loadtxt("Linear_fatigue.out")
	  
	# Linear
	# ax1.plot(defs_linear[:,0], 100*(defs_linear[:,1]-1*defs_linear[0,1]),'b',label = 'Linear',alpha=0.7)

	# GNL
	ax1.plot(defs_gnl[:,0], 100*(defs_gnl[:,1]-1*defs_gnl[0,1]),'r',label = 'GNL',alpha=0.5)


	# start_time = 1
	start_time = 0
	# input_path = '/home/rodrigo/FelipeElgueta/OpenSees-IGA-Examples/modeloAspaFondef/experimental/iga_inputs/case01_190715_13_20_05_y_tip.txt'
	input_path = '/home/rodrigo/FelipeElgueta/OpenSees-IGA-Examples/modeloAspaFondef/experimental/iga_inputs/case00_190609_11_21_47_y_tip.txt'

	dispY_input = np.loadtxt(input_path,skiprows=1)
	nData_input = 4204
	dT_input = 0.00333
	t_input = np.arange(0,nData_input*dT_input,dT_input)
	dispY_input_int = interp1d(t_input,dispY_input)
	dispY_input_interpolated = dispY_input_int(defs_linear[:,0])


	t = defs_gnl[:,0] + start_time
	dispY_input_interpolated = dispY_input_int(t)
	ax1.plot(t-start_time,100*dispY_input_interpolated, '-g', label = 'Measured displacement')

	# # Base displacements
	# ax1.plot(defs_linear[:,0], 100*np.array([0.14/2]*len(defs_linear[:,0])),'--g',label = 'Base displacement')
	# ax1.plot(defs_linear[:,0], 100*np.array([-0.14/2]*len(defs_linear[:,0])),'--g')

	ax1.legend()
	# ax1.ylabel('Y disp [cm]')
	ax1.set_ylabel('Y disp [cm]')
	# ax1.xlabel('Time [s]')
	ax1.set_xlabel('Time [s]')
	ax1.set_title("4.43 Hz fatigue test, vertical displacement at tip")
	ax1.grid()


	# Getting data from input file
	# start_time = 1
	start_time = 0
	# subplot(3,1,2)
	ax2 = fig.add_subplot(3,1,2)

	# input_path = '/home/rodrigo/FelipeElgueta/OpenSees-IGA-Examples/modeloAspaFondef/experimental/iga_inputs/case01_190715_13_20_05_y.txt'
	input_path = '/home/rodrigo/FelipeElgueta/OpenSees-IGA-Examples/modeloAspaFondef/experimental/iga_inputs/case00_190609_11_21_47_y.txt'
	dispY_input = np.loadtxt(input_path,skiprows=1)
	nData_input = 4204
	dT_input = 0.00333
	t_input = np.arange(0,nData_input*dT_input,dT_input)


	dispY_input_int = interp1d(t_input,dispY_input)
	# dispY_input_interpolated = dispY_input_int(defs_linear[:,0])
	dispY_input_interpolated = dispY_input_int(t_input)

	t = defs_gnl[:,0] + start_time

	# t -= start_time
	dispY_input_interpolated = dispY_input_int(t)

	dispYpp=np.gradient(dispY_input_interpolated)


	ax2.plot(t,100*dispY_input_interpolated, label = 'Base displacement')
	ax2.plot(t,100*dispYpp,label = 'acceleration')
	# ax2.ylabel('Base disp [cm]')
	ax2.set_ylabel('Base disp [cm]')
	ax2.set_xlabel('Time [s]')
	ax2.legend()
	ax2.grid()

	# Base displacement history
	# subplot(3,1,3)
	ax3 = fig.add_subplot(3,1,3)

	ax3.plot(t_input,100*dispY_input)
	# t_input -= start_time
	# ax3.plot(t,100*dispY_input_int(t))
	yMin = -100*max(abs(dispY_input))
	ax3.vlines(defs_linear[:,0][-1]+start_time, ymin = yMin, ymax = -yMin, linestyles = 'dashed', colors = 'red', label = 'elapsed time linear')
	ax3.vlines(defs_gnl[:,0][-1]+start_time, ymin = yMin, ymax = -yMin, linestyles = 'dashed', colors = 'green', label = 'elapsed time GNL')

	ax3.set_ylabel('Base disp \n history [cm]')
	ax3.set_xlabel('Time [s]')
	ax3.legend()

	tight_layout()
	ax3.grid()


	# show()


import matplotlib.animation as animation

ani = animation.FuncAnimation(fig,plotAndUpdate,interval=10*1000)
show()