import scipy as sp
import numpy as np

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt

fig = plt.figure(constrained_layout=True)
spec = gridspec.GridSpec(ncols=2, nrows=3, figure=fig)

colors = ['r','b','g','c','m','y']

ω = 4.43 # Hz

# 1, 3, 4 and 6 flapwise

# base case
case_0_data=np.loadtxt('case_0/modalResults_0.txt') # base
case_01_data=np.loadtxt('case_01/modalResults_0.txt') # 1.5 C1
case_02_data=np.loadtxt('case_02/modalResults_0.txt') # 0.5 C1
case_03_data=np.loadtxt('case_03/modalResults_0.txt') # 1.25 C2
case_04_data=np.loadtxt('case_04/modalResults_0.txt') # 0.75 C2
case_05_data=np.loadtxt('case_05/modalResults_0.txt') # 1.25 C3
case_06_data=np.loadtxt('case_06/modalResults_0.txt') # 0.75 C3
case_07_data=np.loadtxt('case_07/modalResults_0.txt') # 1.25 C4
case_08_data=np.loadtxt('case_08/modalResults_0.txt') # 0.75 C4
case_09_data=np.loadtxt('case_09/modalResults_0.txt')
case_10_data=np.loadtxt('case_10/modalResults_0.txt')


# Base case
cycles=case_0_data[:,0]
nModes = len(case_0_data[0,1:])

times=cycles*ω/3600.0/24.0


iColor = 0
iExp = 0
# 1,2,3,4,5,6 modes

# C1 case
ax = fig.add_subplot(spec[0,:])
plt.title("C1 +- 50%")
plt.ylabel("$N_{th}$ freq / $1^{st}$ freq")
plt.xlabel("N cycles")
plt.grid()

from matplotlib.lines import Line2D

line_base = Line2D([0], [0], label='Base Case', color='k', linestyle='-')
line_pos = Line2D([0], [0], label='$+$ C.', color='k', linestyle='--')
line_neg = Line2D([0], [0], label='$-$ C.', color='k', linestyle=':')

plt.legend(handles=[line_base, line_pos, line_neg])

# Base case
freqs = case_0_data[:,1:]
for i in range(nModes):
	if i not in [0,2,3,5]:
		continue
	first_freq=freqs[0,i]

	ax.plot(cycles[:],freqs[:,i]/first_freq, f'-{colors[iColor]}',label=f'Mode {i+1}: {first_freq:.1f} Hz') # Ratio
	# ax.plot(cycles[:],freqs[:,i], f'-{colors[iColor]}x',label=f'Mode {i+1}: {first_freq:.1f} Hz') # Normal
	iColor += 1

# 1.25 C1 case
iColor = 0
cycles=case_01_data[:,0]
nModes = len(case_01_data[0,1:])
freqs = case_01_data[:,1:]
for i in range(nModes):
	if i not in [0,2,3,5]:
		continue
	first_freq=freqs[0,i]	
	ax.plot(cycles[:],freqs[:,i]/first_freq, f'--{colors[iColor]}',label=f'Mode {i+1}: {first_freq:.1f} Hz') # Ratio
	# ax.plot(cycles[:],freqs[:,i], f'-{colors[iColor]}o',label=f'Mode {i+1}: {first_freq:.1f} Hz') # Normal
	iColor += 1

# 0.75 C1 case
iColor = 0
cycles=case_02_data[:,0]
nModes = len(case_02_data[0,1:])
freqs = case_02_data[:,1:]
for i in range(nModes):
	if i not in [0,2,3,5]:
		continue
	first_freq=freqs[0,i]	
	ax.plot(cycles[:],freqs[:,i]/first_freq, f':{colors[iColor]}',label=f'Mode {i+1}: {first_freq:.1f} Hz')
	iColor += 1





# C2 case
ax = fig.add_subplot(spec[1,0])
plt.title("C2 +- 25%")
plt.ylabel("$N_{th}$ freq / $1^{st}$ freq")
plt.xlabel("N cycles")
plt.grid()

# Base case
iColor = 0
freqs = case_0_data[:,1:]
cycles = case_0_data[:,0]
nModes = len(case_0_data[0,1:])
for i in range(nModes):
	if i not in [0,2,3,5]:
		continue
	first_freq=freqs[0,i]	
	ax.plot(cycles[:],freqs[:,i]/first_freq, f'-{colors[iColor]}',label=f'Mode {i+1}: {first_freq:.1f} Hz') # Ratio
	# ax.plot(cycles[:],freqs[:,i], f'-{colors[iColor]}x',label=f'Mode {i+1}: {first_freq:.1f} Hz') # Normal
	iColor += 1

# 1.25 C2 case
iColor = 0
cycles=case_03_data[:,0]
nModes = len(case_03_data[0,1:])
freqs = case_03_data[:,1:]
for i in range(nModes):
	if i not in [0,2,3,5]:
		continue
	first_freq=freqs[0,i]	
	ax.plot(cycles[:],freqs[:,i]/first_freq, f'--{colors[iColor]}',label=f'Mode {i+1}: {first_freq:.1f} Hz') # Ratio
	# ax.plot(cycles[:],freqs[:,i], f'-{colors[iColor]}o',label=f'Mode {i+1}: {first_freq:.1f} Hz') # Normal
	iColor += 1

# 0.75 C2 case
iColor = 0
cycles=case_04_data[:,0]
nModes = len(case_04_data[0,1:])
freqs = case_04_data[:,1:]
for i in range(nModes):
	if i not in [0,2,3,5]:
		continue
	first_freq=freqs[0,i]	
	ax.plot(cycles[:],freqs[:,i]/first_freq, f':{colors[iColor]}',label=f'Mode {i+1}: {first_freq:.1f} Hz')
	iColor += 1


# C3 case
ax = fig.add_subplot(spec[1,1])
plt.title("C3 +- 25%")
plt.ylabel("$N_{th}$ freq / $1^{st}$ freq")
plt.xlabel("N cycles")
plt.grid()

# Base case
iColor = 0
freqs = case_0_data[:,1:]
cycles = case_0_data[:,0]
nModes = len(case_0_data[0,1:])
for i in range(nModes):
	if i not in [0,2,3,5]:
		continue
	first_freq=freqs[0,i]	
	ax.plot(cycles[:],freqs[:,i]/first_freq, f'-{colors[iColor]}',label=f'Mode {i+1}: {first_freq:.1f} Hz') # Ratio
	# ax.plot(cycles[:],freqs[:,i], f'-{colors[iColor]}x',label=f'Mode {i+1}: {first_freq:.1f} Hz') # Normal
	iColor += 1

# 1.25 C3 case
iColor = 0
cycles=case_05_data[:,0]
nModes = len(case_05_data[0,1:])
freqs = case_05_data[:,1:]
for i in range(nModes):
	if i not in [0,2,3,5]:
		continue
	first_freq=freqs[0,i]	
	ax.plot(cycles[:],freqs[:,i]/first_freq, f'--{colors[iColor]}',label=f'Mode {i+1}: {first_freq:.1f} Hz') # Ratio
	# ax.plot(cycles[:],freqs[:,i], f'-{colors[iColor]}o',label=f'Mode {i+1}: {first_freq:.1f} Hz') # Normal
	iColor += 1

# 0.75 C3 case
iColor = 0
cycles=case_06_data[:,0]
nModes = len(case_06_data[0,1:])
freqs = case_06_data[:,1:]
for i in range(nModes):
	if i not in [0,2,3,5]:
		continue
	first_freq=freqs[0,i]	
	ax.plot(cycles[:],freqs[:,i]/first_freq, f':{colors[iColor]}',label=f'Mode {i+1}: {first_freq:.1f} Hz')
	iColor += 1

# # + 50% 
# case_06_50_data=np.loadtxt('case_06/modalResults_1.txt') # 0.5 C3
# # + 75% 
# case_06_75_data=np.loadtxt('case_06/modalResults_2.txt') # 0.75 C3
# # + 100% 
# case_06_1_data=np.loadtxt('case_06/modalResults_3txt') # 1.0 C3
# # + 150% 
# case_06_15_data=np.loadtxt('case_06/modalResults_4.txt') # 1.5 C3

# iColor = 0
# cycles=case_06_50_data[:,0]
# nModes = len(case_06_50_data[0,1:])
# freqs = case_06_50_data[:,1:]
# for i in range(nModes):
# 	if i not in [0,2,3,5]:
# 		continue
# 	first_freq=freqs[0,i]	
# 	ax.plot(cycles[:],freqs[:,i]/first_freq, f'-.{colors[iColor]}',label=f'Mode {i+1}: {first_freq:.1f} Hz')
# 	iColor += 1

# iColor = 0
# cycles=case_06_75_data[:,0]
# nModes = len(case_06_75_data[0,1:])
# freqs = case_06_75_data[:,1:]
# for i in range(nModes):
# 	if i not in [0,2,3,5]:
# 		continue
# 	first_freq=freqs[0,i]	
# 	ax.plot(cycles[:],freqs[:,i]/first_freq, f'-.{colors[iColor]}',label=f'Mode {i+1}: {first_freq:.1f} Hz')
# 	iColor += 1

# iColor = 0
# cycles=case_06_1_data[:,0]
# nModes = len(case_06_1_data[0,1:])
# freqs = case_06_1_data[:,1:]
# for i in range(nModes):
# 	if i not in [0,2,3,5]:
# 		continue
# 	first_freq=freqs[0,i]	
# 	ax.plot(cycles[:],freqs[:,i]/first_freq, f'-.{colors[iColor]}',label=f'Mode {i+1}: {first_freq:.1f} Hz')
# 	iColor += 1

# iColor = 0
# cycles=case_06_15_data[:,0]
# nModes = len(case_06_15_data[0,1:])
# freqs = case_06_15_data[:,1:]
# for i in range(nModes):
# 	if i not in [0,2,3,5]:
# 		continue
# 	first_freq=freqs[0,i]	
# 	ax.plot(cycles[:],freqs[:,i]/first_freq, f'-.{colors[iColor]}',label=f'Mode {i+1}: {first_freq:.1f} Hz')
# 	iColor += 1


# C4 case
ax = fig.add_subplot(spec[2,0])
plt.title("C4 +- 25%")
plt.ylabel("$N_{th}$ freq / $1^{st}$ freq")
plt.xlabel("N cycles")
plt.grid()

# Base case
iColor = 0
freqs = case_0_data[:,1:]
cycles = case_0_data[:,0]
nModes = len(case_0_data[0,1:])
for i in range(nModes):
	if i not in [0,2,3,5]:
		continue
	first_freq=freqs[0,i]	
	ax.plot(cycles[:],freqs[:,i]/first_freq, f'-{colors[iColor]}',label=f'Mode {i+1}: {first_freq:.1f} Hz') # Ratio
	# ax.plot(cycles[:],freqs[:,i], f'-{colors[iColor]}x',label=f'Mode {i+1}: {first_freq:.1f} Hz') # Normal
	iColor += 1

# 1.25 C4 case
iColor = 0
cycles=case_07_data[:,0]
nModes = len(case_07_data[0,1:])
freqs = case_07_data[:,1:]
for i in range(nModes):
	if i not in [0,2,3,5]:
		continue
	first_freq=freqs[0,i]	
	ax.plot(cycles[:],freqs[:,i]/first_freq, f'--{colors[iColor]}',label=f'Mode {i+1}: {first_freq:.1f} Hz') # Ratio
	# ax.plot(cycles[:],freqs[:,i], f'-{colors[iColor]}o',label=f'Mode {i+1}: {first_freq:.1f} Hz') # Normal
	iColor += 1

# 0.75 C4 case
iColor = 0
cycles=case_08_data[:,0]
nModes = len(case_08_data[0,1:])
freqs = case_08_data[:,1:]
for i in range(nModes):
	if i not in [0,2,3,5]:
		continue
	first_freq=freqs[0,i]	
	ax.plot(cycles[:],freqs[:,i]/first_freq, f':{colors[iColor]}',label=f'Mode {i+1}: {first_freq:.1f} Hz')
	iColor += 1

# C5 case
ax = fig.add_subplot(spec[2,1])
plt.title("C5 +- 25%")
plt.ylabel("$N_{th}$ freq / $1^{st}$ freq")
plt.xlabel("N cycles")
plt.grid()

# Base case
iColor = 0
freqs = case_0_data[:,1:]
cycles = case_0_data[:,0]
nModes = len(case_0_data[0,1:])
for i in range(nModes):
	if i not in [0,2,3,5]:
		continue
	first_freq=freqs[0,i]	
	ax.plot(cycles[:],freqs[:,i]/first_freq, f'-{colors[iColor]}',label=f'Mode {i+1}: {first_freq:.1f} Hz') # Ratio
	# ax.plot(cycles[:],freqs[:,i], f'-{colors[iColor]}x',label=f'Mode {i+1}: {first_freq:.1f} Hz') # Normal
	iColor += 1

# 1.25 C5 case
iColor = 0
cycles=case_09_data[:,0]
nModes = len(case_09_data[0,1:])
freqs = case_09_data[:,1:]
for i in range(nModes):
	if i not in [0,2,3,5]:
		continue
	first_freq=freqs[0,i]	
	ax.plot(cycles[:],freqs[:,i]/first_freq, f'--{colors[iColor]}',label=f'Mode {i+1}: {first_freq:.1f} Hz') # Ratio
	# ax.plot(cycles[:],freqs[:,i], f'-{colors[iColor]}o',label=f'Mode {i+1}: {first_freq:.1f} Hz') # Normal
	iColor += 1

# 0.75 C5 case
iColor = 0
cycles=case_10_data[:,0]
nModes = len(case_10_data[0,1:])
freqs = case_10_data[:,1:]
for i in range(nModes):
	if i not in [0,2,3,5]:
		continue
	first_freq=freqs[0,i]	
	ax.plot(cycles[:],freqs[:,i]/first_freq, f':{colors[iColor]}',label=f'Mode {i+1}: {first_freq:.1f} Hz')
	iColor += 1

# plt.tight_layout()
plt.show()