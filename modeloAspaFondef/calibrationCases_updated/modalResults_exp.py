import scipy as sp
import numpy as np

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt



fig = plt.figure(constrained_layout=True, figsize=[10,6])
spec = gridspec.GridSpec(ncols=2, nrows=2, figure=fig)	

colors = ['r','b','g','c','m','y']

ω = 4.43 # Hz

# 1, 3, 4 and 6 flapwise

case_0_data=np.loadtxt('fullCase/modalResults_base.txt') # base
# case_0_data=np.loadtxt('case_01_GNL/modalResults_1.txt') # gnl


# Data from Pate
# First mode
modesExp = [
	np.array([
	[0.0e+0	   , 9.15137133234052],
	[4.5256e+05, 9.1149 ],
	[4.5712e+05, 9.0498 ],
	[8.0155e+05, 9.075  ],
	[8.0585e+05, 9.0717 ],
	[1.1623e+06, 9.0753 ],
	[1.169e+06, 9.0812 ],
	[1.4985e+06, 9.0693 ],
	[2.0432e+06, 9.0792 ],
	[2.3858e+06, 9.0896 ],
	[2.7609e+06, 9.0819 ],
	[3.1166e+06, 9.098  ],
	[3.5754e+06, 9.1031 ],
	[4.2993e+06, 9.1019 ],
	[4.7132e+06, 9.1096 ],
	[5.447e+06, 9.1056 ],
	[5.8837e+06, 9.0951 ],
	[6.1825e+06, 9.0378 ],
	[6.278e+06, 9.0795 ],
	[6.5728e+06, 9.075  ],
	[6.657e+06, 9.0845 ],
	[6.9591e+06, 9.08   ],
	[7.0433e+06, 9.085  ],
	[7.3931e+06, 9.0582 ],
	[8.1245e+06, 9.05   ],
	[8.481e+06, 9.0473 ],
	[8.8597e+06, 9.0278 ],
	[9.2175e+06, 9.0185 ],
	[9.5233e+06, 9.033  ],
	[9.6203e+06, 9.0217 ],
	[9.9304e+06, 9.0185 ],
	[1.0008e+07, 9.0007 ],
	[1.0305e+07, 8.9982 ],
	[1.0765e+07, 8.9921 ],
	[1.2596e+07, 8.9453 ],
	[1.2683e+07, 8.9241 ],
	[1.3064e+07, 8.9194 ],
	[1.3397e+07, 8.9373 ],
	[1.3668e+07, 8.9468 ],
	[1.4071e+07, 8.9157 ],
	[1.4412e+07, 8.9225 ],
	[1.4865e+07, 8.9099 ],
	[1.5525e+07, 8.9069 ],
	[1.6877e+07, 8.8967 ],
	[1.694e+07, 8.8617 ],
	[1.7209e+07, 8.87   ],
	[1.7263e+07, 8.8891 ],
	[1.7509e+07, 8.8711 ],
	[1.9085e+07, 8.9054 ],
	[1.9484e+07, 8.8915 ],
	[1.9563e+07, 8.8907 ],
	[1.9833e+07, 8.8997 ],
	[2.0207e+07, 8.8834 ],
	[2.031e+07, 8.8475 ],
	[2.0606e+07, 8.8551 ],
	[2.0692e+07, 8.8341 ],
	[2.143e+07, 8.8429 ],
	[2.1789e+07, 8.8397 ],
	[2.2133e+07, 8.8098 ],
	[2.2428e+07, 8.7769 ],
	[2.2761e+07, 8.8066 ],
	[2.2834e+07, 8.7901 ],
	[2.3165e+07, 8.7978 ],
	[2.3211e+07, 8.785  ],
	[2.3535e+07, 8.8121 ],
	[3.3639e+07, 8.6816 ],
	[3.3999e+07, 8.6914 ],
	[3.4375e+07, 8.6542 ],
	[3.5211e+07, 8.6107 ],
	[3.5564e+07, 8.611  ],
	[3.5943e+07, 8.6018 ],
	[3.6398e+07, 8.5836 ],
	[3.673e+07, 8.5356 ],
	[3.7123e+07, 8.5743 ],
	[3.7232e+07, 8.572  ],
	[3.7589e+07, 8.549  ],
	[3.7886e+07, 8.5966 ],
	[3.8185e+07, 8.5558 ],
	[3.8628e+07, 8.5221 ],
	[3.8963e+07, 8.5217 ],
	[3.9341e+07, 8.4796 ],
	[3.9418e+07, 8.4643 ],
	[3.9742e+07, 8.483  ],
	[3.9843e+07, 8.4453 ],
	[4.0225e+07, 8.422  ],
	[4.0619e+07, 8.3681 ],
	]),
	np.array([ 
	[0.0e+0	    , 31.5431830530724 ] ,
	[4.5256e+05 , 31.474           ] ,
	[4.5712e+05 , 31.331           ] ,
	[8.0155e+05 , 31.385           ] ,
	[8.0585e+05 , 31.454           ] ,
	[1.1623e+06 , 31.282           ] ,
	[1.169e+06  , 31.26            ] ,
	[1.4985e+06 , 31.168           ] ,
	[2.0432e+06 , 31.017           ] ,
	[2.3858e+06 , 30.894           ] ,
	[2.7609e+06 , 30.845           ] ,
	[3.1166e+06 , 30.854           ] ,
	[3.5754e+06 , 30.883           ] ,
	[4.2993e+06 , 30.913           ] ,
	[4.7132e+06 , 30.939           ] ,
	[5.447e+06  , 30.852           ] ,
	[5.8837e+06 , 30.829           ] ,
	[6.1825e+06 , 30.64            ] ,
	[6.278e+06  , 30.733           ] ,
	[6.5728e+06 , 30.832           ] ,
	[6.657e+06  , 30.913           ] ,
	[6.9591e+06 , 30.798           ] ,
	[7.0433e+06 , 31.081           ] ,
	[7.3931e+06 , 30.931           ] ,
	[8.1245e+06 , 30.873           ] ,
	[8.481e+06  , 30.906           ] ,
	[8.8597e+06 , 30.92            ] ,
	[9.2175e+06 , 30.786           ] ,
	[9.5233e+06 , 30.891           ] ,
	[9.6203e+06 , 30.825           ] ,
	[9.9304e+06 , 30.86            ] ,
	[1.0008e+07 , 30.805           ] ,
	[1.0305e+07 , 30.805           ] ,
	[1.0765e+07 , 30.781           ] ,
	[1.2596e+07 , 30.64            ] ,
	[1.2683e+07 , 30.522           ] ,
	[1.3064e+07 , 30.579           ] ,
	[1.3397e+07 , 30.396           ] ,
	[1.3668e+07 , 30.447           ] ,
	[1.4071e+07 , 30.399           ] ,
	[1.4412e+07 , 30.521           ] ,
	[1.4865e+07 , 30.445           ] ,
	[1.5525e+07 , 30.326           ] ,
	[1.6877e+07 , 30.345           ] ,
	[1.694e+07  , 30.248           ] ,
	[1.7209e+07 , 30.134           ] ,
	[1.7263e+07 , 30.273           ] ,
	[1.7509e+07 , 30.184           ] ,
	[1.9085e+07 , 30.388           ] ,
	[1.9484e+07 , 30.308           ] ,
	[1.9563e+07 , 30.345           ] ,
	[1.9833e+07 , 30.328           ] ,
	[2.0207e+07 , 30.312           ] ,
	[2.031e+07  , 30.167           ] ,
	[2.0606e+07 , 30.229           ] ,
	[2.0692e+07 , 30.161           ] ,
	[2.143e+07  , 30.142           ] ,
	[2.1789e+07 , 30.197           ] ,
	[2.2133e+07 , 30.102           ] ,
	[2.2428e+07 , 29.925           ] ,
	[2.2761e+07 , 30.025           ] ,
	[2.2834e+07 , 30.029           ] ,
	[2.3165e+07 , 30.041           ] ,
	[2.3211e+07 , 29.999           ] ,
	[2.3535e+07 , 30.087           ] ,
	[3.3639e+07 , 30.003           ] ,
	[3.3999e+07 , 30.003           ] ,
	[3.4375e+07 , 29.863           ] ,
	[3.5211e+07 , 29.717           ] ,
	[3.5564e+07 , 29.737           ] ,
	[3.5943e+07 , 29.697           ] ,
	[3.6398e+07 , 29.649           ] ,
	[3.673e+07  , 29.467           ] ,
	[3.7123e+07 , 29.6             ] ,
	[3.7232e+07 , 29.608           ] ,
	[3.7589e+07 , 29.507           ] ,
	[3.7886e+07 , 29.706           ] ,
	[3.8185e+07 , 29.558           ] ,
	[3.8628e+07 , 29.47            ] ,
	[3.8963e+07 , 29.466           ] ,
	[3.9341e+07 , 29.32            ] ,
	[3.9418e+07 , 29.298           ] ,
	[3.9742e+07 , 29.33            ] ,
	[3.9843e+07 , 29.226           ] ,
	[4.0225e+07 , 29.182           ] ,
	[4.0619e+07 , 28.978           ] ,
	[4.1943e+07 , 28.667           ] ,
	] ) ,

	np.array([
	[0.0e+0	    , 74.3945756955987 ] ,
	[4.5256e+05 , 73.532           ] ,
	[4.5712e+05 , 73.18           ] ,
	[8.0155e+05 , 73.264           ] ,
	[8.0585e+05 , 73.241           ] ,
	[1.1623e+06 , 73.249           ] ,
	[1.169e+06  , 73.332           ] ,
	[1.4985e+06 , 73.124           ] ,
	[2.0432e+06 , 73.378           ] ,
	[2.3858e+06 , 73.175           ] ,
	[2.7609e+06 , 73.131           ] ,
	[3.1166e+06 , 73.267           ] ,
	[3.5754e+06 , 73.107           ] ,
	[4.2993e+06 , 73.451           ] ,
	[4.7132e+06 , 73.415           ] ,
	[5.447e+06  , 73.43           ] ,
	[5.8837e+06 , 73.271           ] ,
	[6.1825e+06 , 72.904           ] ,
	[6.278e+06  , 73.159           ] ,
	[6.5728e+06 , 73.221           ] ,
	[6.657e+06  , 73.365           ] ,
	[6.9591e+06 , 73.189           ] ,
	[7.0433e+06 , 73.485           ] ,
	[7.3931e+06 , 73.421           ] ,
	[8.1245e+06 , 73.157           ] ,
	[8.481e+06  , 73.255           ] ,
	[8.8597e+06 , 73.239           ] ,
	[9.2175e+06 , 72.954           ] ,
	[9.5233e+06 , 73.122           ] ,
	[9.6203e+06 , 73.124           ] ,
	[9.9304e+06 , 73.288           ] ,
	[1.0008e+07 , 72.902           ] ,
	[1.0305e+07 , 72.961           ] ,
	[1.0765e+07 , 73.05           ] ,
	[1.2596e+07 , 72.492           ] ,
	[1.2683e+07 , 72.296           ] ,
	[1.3064e+07 , 72.377           ] ,
	[1.3397e+07 , 72.173           ] ,
	[1.3668e+07 , 72.129           ] ,
	[1.4071e+07 , 71.921           ] ,
	[1.4412e+07 , 71.968           ] ,
	[1.4865e+07 , 71.923           ] ,
	[1.5525e+07 , 71.856           ] ,
	[1.6877e+07 , 71.459           ] ,
	[1.694e+07  , 71.4           ] ,
	[1.7209e+07 , 71.307           ] ,
	[1.7263e+07 , 71.619           ] ,
	[1.7509e+07 , 71.29           ] ,
	[1.9085e+07 , 71.831           ] ,
	[1.9484e+07 , 71.693           ] ,
	[1.9563e+07 , 71.763           ] ,
	[1.9833e+07 , 71.735           ] ,
	[2.0207e+07 , 71.678           ] ,
	[2.031e+07  , 71.388           ] ,
	[2.0606e+07 , 71.44           ] ,
	[2.0692e+07 , 71.298           ] ,
	[2.143e+07  , 71.235           ] ,
	[2.1789e+07 , 71.3           ] ,
	[2.2133e+07 , 71.179           ] ,
	[2.2428e+07 , 70.527           ] ,
	[2.2761e+07 , 70.829           ] ,
	[2.2834e+07 , 70.834           ] ,
	[2.3165e+07 , 70.828           ] ,
	[2.3211e+07 , 70.691           ] ,
	[2.3535e+07 , 70.874           ] ,
	[3.3639e+07 , 70.933           ] ,
	[3.3999e+07 , 70.763           ] ,
	[3.4375e+07 , 70.329           ] ,
	[3.5211e+07 , 69.99           ] ,
	[3.5564e+07 , 69.907           ] ,
	[3.5943e+07 , 69.82           ] ,
	[3.6398e+07 , 69.804           ] ,
	[3.673e+07  , 69.532           ] ,
	[3.7123e+07 , 69.774           ] ,
	[3.7232e+07 , 69.699           ] ,
	[3.7589e+07 , 69.402           ] ,
	[3.7886e+07 , 69.866           ] ,
	[3.8185e+07 , 69.584           ] ,
	[3.8628e+07 , 69.288           ] ,
	[3.8963e+07 , 69.27           ] ,
	[3.9341e+07 , 69.053           ] ,
	[3.9418e+07 , 69.004           ] ,
	[3.9742e+07 , 69.01           ] ,
	[3.9843e+07 , 68.813           ] ,
	[4.0225e+07 , 68.776           ] ,
	[4.0619e+07 , 68.479           ] ,
	]),

	np.array([
	[0.0e+0	    , 124.810486652952 ] ,
	# [4.5256e+05 , 125          ] ,
	[4.5256e+05 , 124.810486652952 ] ,
	[4.5712e+05 , 124.4         ] ,
	[8.0155e+05 , 124.34          ] ,
	[8.0585e+05 , 124.51          ] ,
	[1.1623e+06 , 124.43          ] ,
	[1.169e+06  , 124.59          ] ,
	[1.4985e+06 , 124.32          ] ,
	[2.0432e+06 , 123.12          ] ,
	[2.3858e+06 , 123.21          ] ,
	[2.7609e+06 , 123.11          ] ,
	[3.1166e+06 , 123.3          ] ,
	[3.5754e+06 , 123.2          ] ,
	[4.2993e+06 , 123.1          ] ,
	[4.7132e+06 , 123.37          ] ,
	[5.447e+06  , 123.3         ] ,
	[5.8837e+06 , 123.25          ] ,
	[6.1825e+06 , 122.4          ] ,
	[6.278e+06  , 122.94          ] ,
	[6.5728e+06 , 122.93          ] ,
	[6.657e+06  , 123.28          ] ,
	[6.9591e+06 , 122.91          ] ,
	[7.0433e+06 , 123.09          ] ,
	[7.3931e+06 , 123.07          ] ,
	[8.1245e+06 , 122.89          ] ,
	[8.481e+06  , 122.88          ] ,
	[8.8597e+06 , 123.03          ] ,
	[9.2175e+06 , 122.47          ] ,
	[9.5233e+06 , 122.84          ] ,
	[9.6203e+06 , 122.83          ] ,
	[9.9304e+06 , 122.97          ] ,
	[1.0008e+07 , 122.92          ] ,
	[1.0305e+07 , 122.64          ] ,
	[1.0765e+07 , 122.5         ] ,
	[1.2596e+07 , 121.73          ] ,
	[1.2683e+07 , 121.5          ] ,
	[1.3064e+07 , 121.59          ] ,
	[1.3397e+07 , 121.65          ] ,
	[1.3668e+07 , 121.68          ] ,
	[1.4071e+07 , 121.44          ] ,
	[1.4412e+07 , 121.47          ] ,
	[1.4865e+07 , 121.44          ] ,
	[1.5525e+07 , 121.13          ] ,
	[1.6877e+07 , 121.14          ] ,
	[1.694e+07  , 120.61        ] ,
	[1.7209e+07 , 120.6          ] ,
	[1.7263e+07 , 120.87          ] ,
	[1.7509e+07 , 120.72         ] ,
	[1.9085e+07 , 121.34          ] ,
	[1.9484e+07 , 121.33          ] ,
	[1.9563e+07 , 121.31          ] ,
	[1.9833e+07 , 121.49          ] ,
	[2.0207e+07 , 121.28          ] ,
	[2.031e+07  , 120.8          ] ,
	[2.0606e+07 , 120.79         ] ,
	[2.0692e+07 , 120.55          ] ,
	[2.143e+07  , 120.77          ] ,
	[2.1789e+07 , 120.76        ] ,
	[2.2133e+07 , 120.03          ] ,
	[2.2428e+07 , 119.21          ] ,
	[2.2761e+07 , 119.58          ] ,
	[2.2834e+07 , 119.66          ] ,
	[2.3165e+07 , 119.71          ] ,
	[2.3211e+07 , 119.51          ] ,
	[2.3535e+07 , 119.94          ] ,
	[3.3639e+07 , 118.03          ] ,
	[3.3999e+07 , 118.05          ] ,
	[3.4375e+07 , 117.82          ] ,
	[3.5211e+07 , 117.83         ] ,
	[3.5564e+07 , 118.08          ] ,
	[3.5943e+07 , 117.95         ] ,
	[3.6398e+07 , 117.95          ] ,
	[3.673e+07  , 117.37          ] ,
	[3.7123e+07 , 118.07          ] ,
	[3.7232e+07 , 118.2          ] ,
	[3.7589e+07 , 117.83          ] ,
	[3.7886e+07 , 118.34          ] ,
	[3.8185e+07 , 117.82          ] ,
	[3.8628e+07 , 115.53          ] ,
	[3.8963e+07 , 117.51         ] ,
	[3.9341e+07 , 117.46          ] ,
	[3.9418e+07 , 117.37          ] ,
	[3.9742e+07 , 117.7         ] ,
	[3.9843e+07 , 117.1          ] ,
	[4.0225e+07 , 116.83          ] ,
	[4.0619e+07 , 116.55          ] ,
	]),
]

# Base case
cycles=case_0_data[:,0]
nModes = len(case_0_data[0,1:])

iColor = 0
iExp = 0
# 1,2,3,4,5,6 modes

# 1
ax = fig.add_subplot(spec[0,0])
plt.title("First Mode Flap-wise")
plt.ylabel("$N_{th}$ freq / $1^{st}$ freq")
plt.xlabel("N cycles")
plt.grid()

# Base case
freqs = case_0_data[:,1:]
freqs_exp = modesExp[0]
nMode = 0
# if i not in [0,2,3,5]:
# 	continue
first_freq=freqs[0,0]
first_freq_exp=freqs_exp[0,1]

ax.plot(cycles[:],freqs[:,nMode]/first_freq, f'-{colors[iColor]}',label=f'Mode {nMode+1}: {first_freq:.1f} Hz') # Ratio
ax.plot(freqs_exp[:,0],freqs_exp[:,1]/first_freq_exp, f'x{colors[iColor]}',label=f'Experimental Mode {nMode+1}: {first_freq_exp:.1f} Hz') # Ratio
# ax.plot(cycles[:],freqs[:,i], f'-{colors[iColor]}x',label=f'Mode {i+1}: {first_freq:.1f} Hz') # Normal
iColor += 1
plt.legend()

# 2
ax = fig.add_subplot(spec[0,1])
plt.title("Second Mode Flap-wise")
plt.ylabel("$N_{th}$ freq / $1^{st}$ freq")
plt.xlabel("N cycles")
plt.grid()

# Base case
freqs = case_0_data[:,1:]
freqs_exp = modesExp[1]
nMode = 2
# if i not in [0,2,3,5]:
# 	continue
first_freq=freqs[0,nMode]
first_freq_exp=freqs_exp[0,1]

ax.plot(cycles[:],freqs[:,nMode]/first_freq, f'-{colors[iColor]}',label=f'Mode {nMode+1}: {first_freq:.1f} Hz') # Ratio
ax.plot(freqs_exp[:,0],freqs_exp[:,1]/first_freq_exp, f'x{colors[iColor]}',label=f'Experimental Mode {nMode+1}: {first_freq_exp:.1f} Hz') # Ratio
# ax.plot(cycles[:],freqs[:,i], f'-{colors[iColor]}x',label=f'Mode {i+1}: {first_freq:.1f} Hz') # Normal
iColor += 1
plt.legend()

# 3
ax = fig.add_subplot(spec[1,0])
plt.title("Third Mode Flap-wise")
plt.ylabel("$N_{th}$ freq / $1^{st}$ freq")
plt.xlabel("N cycles")
plt.grid()

# Base case
freqs = case_0_data[:,1:]
freqs_exp = modesExp[2]
nMode = 3
# if i not in [0,2,3,5]:
# 	continue
first_freq=freqs[0,nMode]
first_freq_exp=freqs_exp[0,1]

ax.plot(cycles[:],freqs[:,nMode]/first_freq, f'-{colors[iColor]}',label=f'Mode {nMode+1}: {first_freq:.1f} Hz') # Ratio
ax.plot(freqs_exp[:,0],freqs_exp[:,1]/first_freq_exp, f'x{colors[iColor]}',label=f'Experimental Mode {nMode+1}: {first_freq_exp:.1f} Hz') # Ratio
# ax.plot(cycles[:],freqs[:,i], f'-{colors[iColor]}x',label=f'Mode {i+1}: {first_freq:.1f} Hz') # Normal
iColor += 1
plt.legend()

# 4
ax = fig.add_subplot(spec[1,1])
plt.title("Fourth Mode Flap-wise")
plt.ylabel("$N_{th}$ freq / $1^{st}$ freq")
plt.xlabel("N cycles")
plt.grid()

# Base case
freqs = case_0_data[:,1:]
freqs_exp = modesExp[3]
nMode = 5
# if i not in [0,2,3,5]:
# 	continue
first_freq=freqs[0,nMode]
first_freq_exp=freqs_exp[0,1]

ax.plot(cycles[:],freqs[:,nMode]/first_freq, f'-{colors[iColor]}',label=f'Mode {nMode+1}: {first_freq:.1f} Hz') # Ratio
ax.plot(freqs_exp[:,0],freqs_exp[:,1]/first_freq_exp, f'x{colors[iColor]}',label=f'Experimental Mode {nMode+1}: {first_freq_exp:.1f} Hz') # Ratio
# ax.plot(cycles[:],freqs[:,i], f'-{colors[iColor]}x',label=f'Mode {i+1}: {first_freq:.1f} Hz') # Normal
iColor += 1
plt.legend()

plt.show()