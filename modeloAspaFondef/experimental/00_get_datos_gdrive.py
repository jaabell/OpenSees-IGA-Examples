from matplotlib.pylab import *
from subprocess import call
from ver_datos import ver_datos
from get_datos import get_datos
import glob

fname = "190620_16_16_47.txt"

get_datos(fname)

ver_datos("./datos/"+fname,
	save_this_fig=False,
	show_this_fig=True,
	filter_type='sos',
	tmin=5.,
	tmax=10.)

