import scipy as sp
import numpy as np

from matplotlib.pylab import *
filePath_0='/home/rodrigo/FelipeElgueta/OpenSees-IGA-Examples/modeloAspaFondef/4_43HZ/Sensor_Data/20200123-20210615T205444Z-002/20200123/200122_05_36_54.txt'
filePath_1 = '/home/rodrigo/FelipeElgueta/OpenSees-IGA-Examples/modeloAspaFondef/4_43HZ/Sensor_Data/20200123-20210615T205444Z-002/20200123/200122_05_41_55.txt'
filePath_2 = '/home/rodrigo/FelipeElgueta/OpenSees-IGA-Examples/modeloAspaFondef/4_43HZ/Sensor_Data/20200123-20210615T205444Z-002/20200123/200122_05_46_56.txt'

g = 9.807 #m/s2

Data_0=loadtxt(filePath_0,skiprows=38)
t_0=Data_0[:,17]
ac00e_0 = Data_0[:,18] /g
ac00f_0 = Data_0[:,19] /g
Nt_0 = len(ac00e_0)
dt_0 = t_0[1]-t_0[0]
v_0 = np.zeros(Nt_0)
d_0 = np.zeros(Nt_0)
# integrate acceleration to get velocity
v_0[1:] = np.cumsum(ac00e_0[1:] + ac00e_0[0:-1]) * dt_0 / 2  # m/s
# integrate velocity to get displacement
d_0[1:] = np.cumsum(v_0[1:] + v_0[0:-1]) * dt_0 / 2  # m


Data_1=loadtxt(filePath_1,skiprows=38)
t_1=Data_1[:,17]+t_0[-1]
ac00e_1 = Data_1[:,18] /g
ac00f_1 = Data_1[:,19] /g
Nt_1 = len(ac00e_1)
dt_1 = t_1[1]-t_1[0]
v_1 = np.zeros(Nt_1)
d_1 = np.zeros(Nt_1)
# integrate acceleration to get velocity
v_1[1:] = np.cumsum(ac00e_1[1:] + ac00e_1[0:-1]) * dt_1 / 2  # m/s
# integrate velocity to get displacement
d_1[1:] = np.cumsum(v_1[1:] + v_1[0:-1]) * dt_1 / 2  # m

Data_2=loadtxt(filePath_2,skiprows=38)
t_2=Data_2[:,17]+t_1[-1]
ac00e_2 = Data_2[:,18] /g
ac00f_2 = Data_2[:,19] /g
Nt_2 = len(ac00e_2)
dt_2 = t_2[1]-t_2[0]
v_2 = np.zeros(Nt_2)
d_2 = np.zeros(Nt_2)
# integrate acceleration to get velocity
v_2[1:] = np.cumsum(ac00e_2[1:] + ac00e_2[0:-1]) * dt_2 / 2  # m/s
# integrate velocity to get displacement
d_2[1:] = np.cumsum(v_2[1:] + v_2[0:-1]) * dt_2 / 2  # m



plot(t_0,d_0)
plot(t_1,d_1)
plot(t_2,d_2)
show()