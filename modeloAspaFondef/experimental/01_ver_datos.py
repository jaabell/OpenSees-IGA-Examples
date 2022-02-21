from ver_datos import ver_datos
import glob 


files = glob.glob("./datos/*.txt")
Ncols = 3
Nrows = 5

for fname in files:
	ver_datos(fname)
