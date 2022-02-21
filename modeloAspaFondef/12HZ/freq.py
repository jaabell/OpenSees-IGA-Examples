# Info on moving base



from matplotlib.pylab import *
from scipy.interpolate import interp1d

def generateOscillation(amplitude,t_steady,tMax,ω_min, ω_max,nPoints_accel,nPoints_steady):

	deltaT = 0.005/1000

	t=np.arange(0,t_steady,deltaT)

	Fs=1/deltaT
	# Fs=np.linspace(0,1/deltaT,len(t))



	ω_in = np.linspace(ω_min,ω_max,len(t))
	phase_in = np.cumsum(ω_in/Fs)
	y=np.sin(2*np.pi*phase_in)


	t2=np.arange(t_steady+deltaT,tMax+3*deltaT,deltaT)
	y2=np.sin(2*np.pi*ω_max*(t2-0*deltaT))
	t2-=deltaT


	t=np.concatenate([t,t2])
	y=np.concatenate([y,y2])*amplitude

	dispYpp=np.gradient(y)

	# plot(t,y,'-b')
	# plot(t,dispYpp, '-r')



	f = interp1d(t,y,kind = 'cubic')

	# t_int=(np.linspace(0,1,600)**(1/2))*t[-1]
	t_int1=(np.linspace(0,1,nPoints_accel)**(0.75))*t_steady
	t_int2=np.linspace(t_steady,tMax,(tMax-t_steady)*nPoints_steady)

	t_int=np.concatenate([t_int1,t_int2])

	# print(t)
	# print(t_int)



	y_int=f(t_int)
	deltaT_vector = np.diff(t_int)
	deltaT_vector = np.append(deltaT_vector,[deltaT_vector[-1]])
	# print(len(deltaT_vector),len(t_int))
	 
	# for i in range(len(t_int)):
	# 	print(deltaT_vector[i],t_int[i])

	# plot(t_int,y_int,'og')
	# plot(t_int,y_int,'og')

	# show()

	return [deltaT_vector,t_int,y_int]

amplitude = 0.14 / 2 # 14 cm de carrito
t_steady = 2
tMax = 4
ω_min = 0   # Hz 
ω_max = 10   # Hz
nPoints_accel=400
nPoints_steady=400
[deltaT_vector,t,y] = generateOscillation(amplitude,t_steady,tMax,ω_min, ω_max,nPoints_accel,nPoints_steady)
dispYpp=np.gradient(y)
plot(t,dispYpp, '-r')
plot(t,y,'-b')
print(len(t))
plot(t,y,'og')
show()
print(deltaT_vector)
print(t)