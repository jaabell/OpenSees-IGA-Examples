import scipy as sp
import numpy as np

from matplotlib.pylab import *

data=loadtxt('modalResults_0.txt')

cycles=data[:,0]
nModes = len(data[0,1:])

fig,ax=subplots()
# ax2 = ax.twiny()

# 1, 3, 4 and 6 flapwise

ω = 4.43 #Hz
times=cycles*ω/3600.0/24.0

colors = ['r','b','g','c','m','y']
iColor = 0
iExp = 0
# 1,2,3,4,5,6 modes



freqs = data[:,1:]
for i in range(nModes):

	if i not in [0,2,3,5]:
		continue
		
	first_freq=freqs[0,i]	
	ax.plot(cycles[:],freqs[:,i]/first_freq, f'-{colors[iColor]}x',label=f'Mode {i+1}: {first_freq:.1f} Hz')
	
	iColor += 1



ax.set_xlabel("Cycles")
# ax2.set_xlabel("Days")
ax.set_ylabel("Ratio of frequency")
ax.grid()
ax.legend()
title("Fatigue modal degradation",loc='center',fontsize=20)
tight_layout()
show()