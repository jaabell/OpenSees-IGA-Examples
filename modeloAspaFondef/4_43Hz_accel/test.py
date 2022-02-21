import numpy as np
import scipy as sp

from matplotlib.pylab import *

accel_path =  '/home/rodrigo/FelipeElgueta/OpenSees-IGA-Examples/modeloAspaFondef/experimental/iga_inputs/case00_190609_11_21_47_ypp.txt'
accel_path_unfiltered =  '/home/rodrigo/FelipeElgueta/OpenSees-IGA-Examples/modeloAspaFondef/experimental/iga_inputs/case00_190609_11_21_47_ypp_unfiltered.txt'

ypp = np.loadtxt(accel_path,skiprows=1)
nData_input = len(ypp)
dT_input = 0.00333
t_input = np.arange(0,nData_input*dT_input,dT_input)
print(len(t_input))



from scipy.interpolate import interp1d

ypp_interpolator = interp1d(t_input,ypp)
# dT = 0.005
Fs = 250.0
# dT = 0.01
dT=1.0/Fs
print("dT = ", dT)
# start_time = 1
start_time = 0
t = np.arange(start_time,t_input[-1],dT)
ypp_interpolated = ypp_interpolator(t)

plot(t_input,ypp,'b',label = 'filtered')
plot(t,ypp_interpolated,'og')
# show()



# Unfiltered

ypp = np.loadtxt(accel_path_unfiltered,skiprows=1)
nData_input = len(ypp)
dT_input = 0.00333
t_input = np.arange(0,nData_input*dT_input,dT_input)
print(len(t_input))

plot(t_input,ypp,'r', label = 'original')


legend()
show()
# # FFT

# from scipy.fft import rfft, rfftfreq

# normalized_data = np.int16((ypp / ypp.max()) * 32767)

# yf = rfft(normalized_data)
# xf = rfftfreq(nData_input, 1 / dT_input)

# plt.plot(xf, np.abs(yf))
# plt.show()
